#!/usr/bin/env python3
"""
ORCA ECD Data Processing Script
ORCA ECD数据处理脚本
"""

import os
import re
import numpy as np

# Try importing matplotlib; if failed, only generate data files
# 尝试导入matplotlib，如果失败则只生成数据文件
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not installed, will only generate data files")
    print("警告: matplotlib未安装，将只生成数据文件")

def read_weights(weights_file):
    """Read weight file with multiple encoding attempts
    读取权重文件，尝试多种编码"""
    weights = {}
    
    if not os.path.exists(weights_file):
        print(f"Error: Weight file does not exist: {weights_file}")
        print(f"错误: 权重文件不存在: {weights_file}")
        return weights
    
    # Try multiple encoding methods
    # 尝试多种编码方式
    encodings = ['utf-8', 'gbk', 'gb2312', 'latin-1']
    
    for encoding in encodings:
        try:
            with open(weights_file, 'r', encoding=encoding, errors='ignore') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            conf_id = int(parts[0].strip())
                            weight = float(parts[1].strip())
                            weights[conf_id] = weight
                        except ValueError:
                            continue
            print(f"Successfully read weight file using {encoding} encoding")
            print(f"使用 {encoding} 编码成功读取权重文件")
            break
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(f"Error reading with {encoding} encoding: {e}")
            print(f"使用 {encoding} 编码读取时出错: {e}")
            continue
    
    return weights

def parse_ecd_file(filepath, debug=False):
    """Robust parsing of ORCA ECD output (compatible with multiple column formats), returns (energies_array, R_array) or None
    鲁棒解析 ORCA ECD 输出（兼容多种列格式），返回 (energies_array, R_array) 或 None"""
    import math

    if not os.path.exists(filepath):
        if debug: 
            print(f"parse_ecd_file: File does not exist: {filepath}")
            print(f"parse_ecd_file: 文件不存在: {filepath}")
        return None

    encodings = ['utf-8', 'gbk', 'gb2312', 'latin-1', 'utf-16']
    content = None
    used_encoding = None
    for enc in encodings:
        try:
            with open(filepath, 'r', encoding=enc, errors='ignore') as f:
                content = f.read()
            used_encoding = enc
            break
        except Exception:
            continue

    if content is None:
        print(f"Warning: Unable to read {filepath} with common encodings")
        print(f"警告: 无法以常见编码读取 {filepath}")
        return None

    # Find text range of ECD/CD blocks
    # 查找 ECD/CD 块的文本范围
    # Try several header identifiers
    # 尝试几种 header 标识
    start_patterns = [
        r'CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS', 
        r'CD SPECTRUM',
        r'CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS'
    ]
    start_idx = None
    for p in start_patterns:
        m = re.search(p, content)
        if m:
            start_idx = m.start()
            break

    if start_idx is None:
        if debug:
            print(f"parse_ecd_file: CD SPECTRUM marker not found in {os.path.basename(filepath)} (tried encoding: {used_encoding})")
            print(f"parse_ecd_file: 在 {os.path.basename(filepath)} 中找不到 CD SPECTRUM 标记 (尝试的编码: {used_encoding})")
            print("=== DEBUG START (file head 2000 chars) ===")
            print(content[:2000])
            print("=== DEBUG END ===")
        return None

    # Find block end from start_idx backwards (next '---' or empty line, or 'Total run time', etc.)
    # 从 start_idx 向后找到块结束（下一组 '---' 或空行，或 'Total run time' 等）
    tail_markers = ['CD SPECTRUM VIA TRANSITION VELOCITY', 'Total run time', '\n\n', '\n\n\n']
    end_idx = None
    # prefer explicit separator blocks bounded by --- lines
    sep_match = re.search(r'\n-+\n', content[start_idx:])
    if sep_match:
        # take region between first separator after header and next big separator or next property header
        # We'll take up to next two '---' blocks to be safe
        sec_start = start_idx
        # search for next major header after sec_start+len
        # simple approach: take next 5000 chars chunk starting at start_idx
        chunk = content[start_idx:start_idx+10000]
        # cut at next "CD SPECTRUM VIA" or "Total run time" or "ABSORPTION SPECTRUM" etc.
        for marker in ['CD SPECTRUM VIA TRANSITION VELOCITY','Total run time','ABSORPTION SPECTRUM','PROPERTY CALCULATIONS']:
            mi = chunk.find(marker)
            if mi != -1:
                end_idx = start_idx + mi
                break
        if end_idx is None:
            end_idx = start_idx + len(chunk)
    else:
        end_idx = min(start_idx + 10000, len(content))

    block = content[start_idx:end_idx]

    # Filter lines containing "->" (usually actual data lines)
    # 按行筛出包含"->"的过渡行（一般为实际数据行）
    lines = block.splitlines()
    data_lines = []
    for ln in lines:
        if '->' in ln and re.search(r'\d', ln):
            # Exclude header lines (containing 'Energy', 'Wavelength', etc.) and short lines
            # 排除表头线（含 'Energy', 'Wavelength' 等）和短行
            if len(ln.strip()) < 10:
                continue
            # ignore lines that are clearly separators
            if re.match(r'^\s*-+\s*$', ln):
                continue
            data_lines.append(ln)

    if not data_lines:
        if debug:
            print(f"parse_ecd_file: No data lines found in block (file {os.path.basename(filepath)})")
            print(f"parse_ecd_file: 未在块中发现数据行 (文件 {os.path.basename(filepath)})")
            print("=== DEBUG BLOCK START ===")
            print(block[:2000])
            print("=== DEBUG BLOCK END ===")
        return None

    energies = []
    R_values = []

    for ln in data_lines:
        # Extract all tokens that can be converted to float
        # 提取所有可以转换为 float 的 token
        toks = ln.split()
        floats = []
        for t in toks:
            try:
                # Filter out arrows and letter combinations
                # 过滤掉箭头和字母组合
                if re.match(r'^-?[\d]+\.[\d]+$', t) or re.match(r'^-?\d+\.\d+[eE][\+\-]?\d+$', t):
                    vf = float(t)
                    if math.isfinite(vf):
                        floats.append(vf)
                else:
                    # try to strip trailing punctuation
                    tt = t.strip().strip(',')
                    vf = float(tt)
                    floats.append(vf)
            except Exception:
                continue

        # If no floats detected, try extracting substrings with decimals from line
        # 如果没有检测到 floats，则尝试从 line 中提取带小数的子串
        if not floats:
            found = re.findall(r'-?\d+\.\d+(?:[eE][\+\-]?\d+)?', ln)
            floats = [float(x) for x in found]

        if len(floats) < 2:
            # Invalid line, skip
            # 无效行，跳过
            if debug:
                print(f"parse_ecd_file: Skipping line (insufficient floats): {ln}")
                print(f"parse_ecd_file: 跳过行（float不足）: {ln}")
            continue

        # Usually floats list is:
        # [energy(eV), wavenumber(cm-1), wavelength(nm), R, MX, MY, MZ]
        # Common index: energy = floats[0]
        # Most robust strategy: energy = first float in reasonable eV range (0.5 - 30)
        # 通常 floats 列表为:
        # [energy(eV), wavenumber(cm-1), wavelength(nm), R, MX, MY, MZ]
        # 常见索引: energy = floats[0]
        # 最稳妥的策略：energy = first float in reasonable eV range (0.5 - 30)
        energy = None
        for f in floats:
            if 0.1 < abs(f) < 30:  # energy in eV
                energy = float(f)
                break
        if energy is None:
            energy = floats[0]

        # R position is usually the 4th float (index 3) in ORCA outputs above; if not present, take first negative/large value
        # R 的位置通常是 the 4th float (index 3) in ORCA outputs above; 若不存在则取第一个负/较大值
        R = None
        if len(floats) >= 4:
            R = floats[3]
        else:
            # Fallback strategy: find first float with large absolute value (excluding energy/wavenumber/wavelength)
            # 回退策略：找第一个绝对值较大的 float（排除能量/波数/波长）
            cand = [f for f in floats if abs(f) > 1e-3 and (f < -0.5 or f > 0.5)]
            if cand:
                R = cand[0]
            else:
                # Finally, take item near the end floats[-(number)]
                # 最后退化取 floats[-(number)] 末尾靠近的项
                R = floats[-4] if len(floats) >= 4 else floats[-1]

        # Validate energy and R seem reasonable (energy positive, R can be positive or negative but not NaN)
        # 验证 energy 与 R 看起来合理（energy 正，R 可正可负但不为 NaN）
        try:
            energy = float(energy)
            R = float(R)
            energies.append(energy)
            R_values.append(R)
        except Exception:
            if debug:
                print(f"parse_ecd_file: Conversion failed ln: {ln} floats: {floats}")
                print(f"parse_ecd_file: 转换失败 ln: {ln} floats: {floats}")

    if not energies:
        if debug:
            print(f"parse_ecd_file: No valid energies after parsing (file {os.path.basename(filepath)})")
            print(f"parse_ecd_file: 解析后无有效能量 (file {os.path.basename(filepath)})")
        return None

    if debug:
        print(f"parse_ecd_file: Parsed {len(energies)} transitions in {os.path.basename(filepath)} (encoding={used_encoding})")
        print(f"parse_ecd_file: 在 {os.path.basename(filepath)} 中解析到 {len(energies)} 个跃迁 (encoding={used_encoding})")

    return np.array(energies), np.array(R_values)

def calculate_weighted_spectrum(opt_dir, ecd_dir, sigma=0.3):
    """Calculate weighted ECD spectrum
    计算加权ECD谱"""
    print("Starting calculation of weighted ECD spectrum...")
    print("开始计算加权ECD谱...")
    
    # Try reading weight file
    # 尝试读取权重文件
    weights_file = os.path.join(opt_dir, 'ecd_weights.txt')
    weights = read_weights(weights_file)
    
    if not weights:
        print("Error: Unable to obtain weight data")
        print("错误: 无法获取权重数据")
        return None, None
    
    print(f"Read weights for {len(weights)} conformations")
    print(f"读取到 {len(weights)} 个构象的权重")
    
    # Collect all data
    # 收集所有数据
    all_data = {}
    
    for conf_id, weight in weights.items():
        filename = f"ecd-conf-{conf_id}.out"
        filepath = os.path.join(ecd_dir, filename)
        
        data = parse_ecd_file(filepath)
        if data is not None:
            all_data[conf_id] = {
                'energy': data[0],
                'R': data[1],
                'weight': weight
            }
        else:
            print(f"Warning: Unable to parse data for conformation {conf_id}")
            print(f"警告: 无法解析构象 {conf_id} 的数据")
    
    if not all_data:
        print("Error: No available ECD data")
        print("错误: 没有可用的ECD数据")
        return None, None
    
    print(f"Successfully parsed ECD data for {len(all_data)} conformations")
    print(f"成功解析 {len(all_data)} 个构象的ECD数据")
    
    # Create energy grid (1.5-8.5 eV, corresponding to approximately 146-827 nm)
    # 创建能量网格 (1.5-8.5 eV，对应约146-827 nm)
    energy_grid = np.linspace(1.5, 8.5, 1000)
    spectrum = np.zeros_like(energy_grid)
    
    # Calculate weighted average spectrum
    # 计算加权平均谱
    total_weight = sum(item['weight'] for item in all_data.values())
    
    for conf_id, conf_data in all_data.items():
        weight = conf_data['weight'] / total_weight
        
        for energy, R in zip(conf_data['energy'], conf_data['R']):
            # Gaussian broadening
            # 高斯展宽
            gaussian = np.exp(-(energy_grid - energy)**2 / (2 * sigma**2))
            gaussian /= np.sqrt(2 * np.pi * sigma**2)
            spectrum += weight * R * gaussian
    
    return energy_grid, spectrum

def fft_smooth(y, smooth_factor):
    """
    Smooth data using FFT low-pass filtering
    使用 FFT 低通滤波平滑数据
    smooth_factor: 0~1, smaller is smoother, recommended 0.02~0.2
    smooth_factor: 0~1，越小越平滑，推荐 0.02~0.2
    """
    y = np.array(y)
    n = len(y)

    # FFT
    Y = np.fft.rfft(y)

    # Frequency cutoff (low-pass)
    # 频率裁剪（低通）
    cutoff = int(len(Y) * smooth_factor)   # Keep first smooth_factor percent of frequencies
    Y[cutoff:] = 0                         # Set high frequencies to zero
    
    # Inverse FFT
    # 逆 FFT
    y_smooth = np.fft.irfft(Y, n=n)
    return y_smooth

def read_experimental_data(exp_file='result.csv'):
    """Read experimental data file (supports two formats: with XYDATA marker and simple CSV format)
    读取实验数据文件（支持两种格式：带XYDATA标记的格式和简单CSV格式）"""
    if not os.path.exists(exp_file):
        print(f"Warning: Experimental data file does not exist: {exp_file}")
        print(f"警告: 实验数据文件不存在: {exp_file}")
        return None, None
    
    wavelengths = []
    intensities = []
    
    # Try multiple encodings
    # 尝试多种编码
    encodings = ['utf-8', 'gbk', 'gb2312', 'latin-1']
    
    for encoding in encodings:
        try:
            with open(exp_file, 'r', encoding=encoding) as f:
                lines = f.readlines()
            
            # Check if file contains XYDATA marker
            # 检查文件是否包含XYDATA标记
            has_xydata = any('XYDATA' in line for line in lines)
            
            if has_xydata:
                # Format 1: Format with XYDATA marker
                # 格式1: 包含XYDATA标记的格式
                # Find data start position
                # 寻找数据开始的位置
                start_idx = -1
                for i, line in enumerate(lines):
                    if 'XYDATA' in line:
                        start_idx = i + 1
                        break
                
                if start_idx == -1:
                    print("Warning: XYDATA marker not found in experimental data file")
                    print("警告: 实验数据文件中未找到XYDATA标记")
                    return None, None
                
                # Read data
                # 读取数据
                for line in lines[start_idx:]:
                    # Skip empty lines and comment lines
                    # 跳过空行和注释行
                    if not line.strip() or ',' not in line:
                        continue
                    
                    # Check if reached comment section
                    # 检查是否到达注释部分
                    if '##### Extended Information' in line or '[Comments]' in line:
                        break
                    
                    # Parse data line
                    # 解析数据行
                    parts = line.strip().split(',')
                    if len(parts) >= 2:
                        try:
                            wavelength = float(parts[0].strip())
                            intensity = float(parts[1].strip())
                            wavelengths.append(wavelength)
                            intensities.append(intensity)
                        except ValueError:
                            continue
            else:
                # Format 2: Simple CSV format (each line: wavelength,intensity)
                # 格式2: 简单CSV格式（每行: 波长,强度）
                for line in lines:
                    # Skip empty lines and comment lines (starting with #)
                    # 跳过空行和注释行（以#开头的行）
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Try to parse data line
                    # 尝试解析数据行
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            # Handle scientific notation (e.g., 2.18646E-05)
                            # 处理科学计数法（如2.18646E-05）
                            wavelength_str = parts[0].strip()
                            intensity_str = parts[1].strip()
                            
                            # Replace E/e in scientific notation with Python-recognizable format
                            # 替换科学计数法中的E/e为python可识别的格式
                            if 'E' in intensity_str or 'e' in intensity_str:
                                # Use float() directly, it can handle scientific notation
                                # 直接使用float()，它能处理科学计数法
                                intensity = float(intensity_str)
                            else:
                                intensity = float(intensity_str)
                            
                            wavelength = float(wavelength_str)
                            
                            wavelengths.append(wavelength)
                            intensities.append(intensity)
                        except ValueError as e:
                            print(f"Error parsing line '{line}': {e}")
                            print(f"解析行时出错 '{line}': {e}")
                            continue
            
            if wavelengths:  # If data read successfully
                print(f"Successfully read experimental data file using {encoding} encoding")
                print(f"使用 {encoding} 编码成功读取实验数据文件")
                break
                
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(f"Error reading experimental data with {encoding} encoding: {e}")
            print(f"使用 {encoding} 编码读取实验数据时出错: {e}")
            continue
    
    if not wavelengths:
        print("Warning: No experimental data read")
        print("警告: 未读取到实验数据")
        return None, None
    
    # Convert to numpy arrays and ensure correct order (from low to high wavelength)
    # 转换为numpy数组并确保顺序正确（从低波长到高波长）
    wavelengths = np.array(wavelengths)
    intensities = np.array(intensities)
    
    # Check data order: if wavelengths are descending (e.g., 400, 399.9, 399.8...), reverse
    # 检查数据顺序：如果波长是降序（如400, 399.9, 399.8...），则反转
    if len(wavelengths) > 1 and wavelengths[0] > wavelengths[-1]:
        print("Detected descending wavelength order, reversing data...")
        print("检测到波长降序排列，正在反转数据...")
        wavelengths = wavelengths[::-1]
        intensities = intensities[::-1]
    
    print(f"Read {len(wavelengths)} experimental data points")
    print(f"读取到 {len(wavelengths)} 个实验数据点")
    print(f"Wavelength range: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} nm")
    print(f"波长范围: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} nm")
    print(f"Intensity range: {intensities.min():.6f} - {intensities.max():.6f}")
    print(f"强度范围: {intensities.min():.6f} - {intensities.max():.6f}")
    
    return wavelengths, intensities

def save_combined_spectrum_csv(calc_wavelength_grid, calc_spectrum, exp_wavelengths, exp_intensities, sigma, shift, scale_factor, smooth_factor=None):
    """Save combined spectrum data to CSV, including calculated data (normalized and scaled) and experimental data (after smoothing)
    保存合并的谱图数据到CSV，包含计算数据（归一化并缩放）和实验数据（平滑后）"""
    filename = f'ecd_combined_spectrum_sigma{sigma}_shift{shift}_scale{scale_factor}_smooth{smooth_factor}.csv'
    print(f"Saving combined CSV data to {filename}...")
    print(f"保存合并的CSV数据到 {filename}...")
    
    # ===== Normalization of calculated data =====
    # ===== 计算数据的归一化 =====
    max_calc_intensity = np.max(np.abs(calc_spectrum))
    if max_calc_intensity == 0:
        norm_calc_spectrum = calc_spectrum
    else:
        norm_calc_spectrum = calc_spectrum / max_calc_intensity
    
    # ===== Apply scaling factor =====
    # ===== 应用缩放因子 =====
    scaled_calc_spectrum = norm_calc_spectrum * scale_factor
    
    # ===== Generate inverted data =====
    # ===== 生成倒置数据 =====
    inverted_calc_spectrum = -scaled_calc_spectrum
    # ======================
    
    # Write CSV file
    # 写入CSV文件
    with open(filename, 'w', encoding='utf-8') as f:
        # Write file header
        # 写入文件头
        f.write("# ECD Spectrum Data - Calculated data (normalized and scaled) vs Experimental data (smoothed)\n")
        f.write("# ECD谱数据 - 计算数据（归一化并缩放）与实验数据（平滑）对比\n")
        f.write("# Note: Calculated data normalized (max intensity = 1) and multiplied by scaling factor{}\n".format(scale_factor))
        f.write("# 注意：计算数据已归一化（最大强度 = 1）并乘以缩放因子{}\n".format(scale_factor))
        f.write("# Smoothed experimental data (mdeg units)\n")
        f.write("# 平滑后实验数据（mdeg单位）\n")
        f.write("# sigma = {} eV, shift = {} eV, scale_factor = {}, smooth_factor = {}\n".format(sigma, shift, scale_factor, smooth_factor))
        f.write("\n")
        
        # Write calculated data title
        # 写入计算数据标题
        f.write("=== Calculated Data (normalized, scaled, and inverted version) ===\n")
        f.write("=== 计算数据（归一化并缩放，及倒置版本）===\n")
        f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted\n")
        
        # Write calculated data
        # 写入计算数据
        for wl, scaled_inten, inverted_inten in zip(calc_wavelength_grid, scaled_calc_spectrum, inverted_calc_spectrum):
            f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f}\n")
        
        f.write("\n")
        
        # Write experimental data title
        # 写入实验数据标题
        f.write("=== Experimental Data (smoothed) ===\n")
        f.write("=== 实验数据（平滑）===\n")
        f.write("Wavelength(nm),Experimental_Intensity(mdeg)\n")
        
        # Write experimental data
        # 写入实验数据
        for wl, inten in zip(exp_wavelengths, exp_intensities):
            f.write(f"{wl:.4f},{inten:.6f}\n")
    
    print(f"Combined CSV data saved: {filename}")
    print(f"合并的CSV数据保存完成: {filename}")
    
    # Also save an Origin-friendly format (aligned columns)
    # 同时保存一个Origin友好格式（对齐的列）
    origin_filename = f'ecd_origin_format_sigma{sigma}_shift{shift}_scale{scale_factor}_smooth{smooth_factor}.csv'
    print(f"Saving Origin-format CSV data to {origin_filename}...")
    print(f"保存Origin格式CSV数据到 {origin_filename}...")
    
    # Create interpolation function for experimental data (linear interpolation)
    # 创建实验数据的插值函数（线性插值）
    try:
        # Check if there are enough experimental data points for interpolation
        # 检查是否有足够的实验数据点进行插值
        if len(exp_wavelengths) > 1:
            # Manually implement linear interpolation (avoid scipy dependency)
            # 手动实现线性插值（避免依赖scipy）
            exp_interp_intensities = np.interp(calc_wavelength_grid, exp_wavelengths, exp_intensities, 
                                               left=np.nan, right=np.nan)
            
            # Write Origin-friendly format
            # 写入Origin友好格式
            with open(origin_filename, 'w', encoding='utf-8') as f:
                f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted,Experimental_Intensity(mdeg)\n")
                for i, wl in enumerate(calc_wavelength_grid):
                    scaled_inten = scaled_calc_spectrum[i]
                    inverted_inten = inverted_calc_spectrum[i]
                    exp_inten = exp_interp_intensities[i]
                    f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f},{exp_inten:.6f}\n")
            
            print(f"Origin-format CSV data saved: {origin_filename}")
            print(f"Origin格式CSV数据保存完成: {origin_filename}")
        else:
            print("Too few experimental data points, cannot interpolate, skipping Origin format save")
            print("实验数据点太少，无法进行插值，跳过Origin格式保存")
    except Exception as e:
        print(f"Interpolation failed, skipping Origin format save: {e}")
        print(f"插值失败，跳过Origin格式保存: {e}")
        
        # Save simplified version
        # 保存简化版本
        with open(origin_filename, 'w', encoding='utf-8') as f:
            f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted\n")
            for i, wl in enumerate(calc_wavelength_grid):
                scaled_inten = scaled_calc_spectrum[i]
                inverted_inten = inverted_calc_spectrum[i]
                f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f}\n")
        print("Saved simplified version with only calculated data")
        print("已保存仅含计算数据的简化版本")

def plot_spectrum(energy_grid, calc_spectrum, sigma, shift, scale_factor, exp_wavelengths=None, exp_intensities=None, plot_exp_name=None, plot_calc_name=None, plot_calc_name_invert=None, smooth_factor=None):
    if not HAS_MATPLOTLIB:
        print("matplotlib not installed, skipping plotting")
        print("matplotlib未安装，跳过绘图")
        return
    
    # Set Times New Roman font
    # 设置Times New Roman字体
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'  # Math font also set to Times New Roman style
    
    # ===== Normalization and scaling of calculated data =====
    # ===== 计算数据的归一化和缩放 =====
    max_calc_intensity = np.max(np.abs(calc_spectrum))
    if max_calc_intensity == 0:
        norm_calc_spectrum = calc_spectrum
    else:
        norm_calc_spectrum = calc_spectrum / max_calc_intensity
    
    # Apply scaling factor
    # 应用缩放因子
    scaled_calc_spectrum = norm_calc_spectrum * scale_factor
    
    # Generate inverted data
    # 生成倒置数据
    inverted_calc_spectrum = -scaled_calc_spectrum
    # ======================
    
    plt.figure(figsize=(10, 6))
    
    # Calculated data: main coordinate changed to wavelength
    # 计算数据：主坐标改为波长
    wavelength_grid = 1239.84 / energy_grid
    
    # Plot calculated spectrum (scaled)
    # 绘制计算谱图（缩放后）
    plt.plot(wavelength_grid, scaled_calc_spectrum, 'r--', linewidth=1.5, 
             label=f'{plot_calc_name}')
    
    # Plot inverted calculated spectrum
    # 绘制倒置计算谱图
    plt.plot(wavelength_grid, inverted_calc_spectrum, 'b--', linewidth=1.5, alpha=0.7,
             label=f'{plot_calc_name_invert}')
    
    # Plot experimental data (if provided)
    # 绘制实验数据（如果提供）
    if exp_wavelengths is not None and exp_intensities is not None:
        # Experimental data not normalized, keep original values
        # 实验数据不进行归一化，保持原始值
        plt.plot(exp_wavelengths, exp_intensities, 'k-', linewidth=1.5, label=f'{plot_exp_name}')
    
    # Main x-axis: 200–400 nm, uniform scale
    # 主 x 轴：200–400 nm，均匀刻度
    plt.xlim(200, 400)
    plt.xticks([200, 250, 300, 350, 400])
    plt.xlabel('Wavelength (nm)', fontsize=14, fontweight='bold')
    
    # y-axis
    # y 轴
    plt.ylabel('Δε (M$^{-1}$cm$^{-1}$)', fontsize=14, fontweight='bold')
    
    # Set title based on whether experimental data is present
    # 根据是否有实验数据设置标题
    if exp_wavelengths is not None and exp_intensities is not None:
        plt.title(f'ECD Spectrum Comparison\nσ={sigma} eV, shift={shift} eV, scale={scale_factor}, smooth={smooth_factor}', fontsize=14)
    else:
        plt.title(f'Calculated ECD Spectrum\nσ={sigma} eV, shift={shift} eV, scale={scale_factor}', fontsize=14)
    
    # Zero line
    # 零线
    plt.axhline(y=0, color='black', linewidth=0.5)
    
    # Secondary x-axis (top): energy eV
    # 副 x 轴（上方）：能量 eV
    ax = plt.gca()
    
    def nm_to_ev(x):
        return 1239.84 / x
    
    def ev_to_nm(x):
        return 1239.84 / x
    
    ax2 = ax.secondary_xaxis('top', functions=(nm_to_ev, ev_to_nm))
    ax2.set_xlabel('Energy (eV)', fontsize=12)
    
    # Add legend
    # 添加图例
    plt.legend(loc='best', shadow=True, edgecolor='black', facecolor='white', fontsize=10)
    
    # Grid
    # 网格
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Set output filename based on whether experimental data is present
    # 根据是否有实验数据设置输出文件名
    if exp_wavelengths is not None and exp_intensities is not None:
        output_file = f'ecd_spectrum_comparison_sigma{sigma}_shift{shift}_scale{scale_factor}_smooth{smooth_factor}.png'
    else:
        output_file = f'ecd_spectrum_sigma{sigma}_shift{shift}_scale{scale_factor}.png'
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Spectrum saved as {output_file}")
    print(f"谱图已保存为 {output_file}")
    plt.show()

def main():
    """Main function
    主函数"""
    print("="*60)
    print("ORCA ECD Data Processing Script")
    print("ORCA ECD数据处理脚本")
    print("="*60)
    
    # Set calculation parameters
    # 设置计算参数
    opt_dir = 'opt_conf'    # Optimization results folder
    ecd_dir = 'ecd_opt'     # ECD results folder
    sigma = 0.1             # Gaussian broadening parameter
    shift = 0.0             # Spectrum shift (negative shifts right, e.g., -0.2)
    scale_factor = 1        # Scaling factor for calculated data (for alignment with experimental data)

    # Set experimental parameters
    # 设置实验参数
    exp_file = 'exp_data.csv' # Experimental data file
    smooth_factor = 0.05      # FFT factor for smoothing experimental data (0~1, smaller is smoother)

    # Set legend names
    # 设置图例名称
    plot_exp_name = 'Experimental ECD of your_molecule_name'
    plot_calc_name = "Calculated ECD of ($R$)- your_molecule_name"
    plot_calc_name_invert = "Calculated ECD of ($S$)- your_molecule_name"
    
    # Get parameters from user
    # 从用户获取参数
    print("Current parameter settings:")
    print("当前参数设置：")
    print(f"  sigma (Gaussian broadening): {sigma}")
    print(f"  sigma (高斯展宽): {sigma}")
    print(f"  shift (spectrum shift): {shift}")
    print(f"  shift (能谱平移): {shift}")
    print(f"  scale_factor (calculated data scaling factor): {scale_factor}")
    print(f"  scale_factor (计算数据缩放因子): {scale_factor}")
    print(f"  Smoothing factor (FFT smoothing for experimental data): {smooth_factor}")
    print(f"  平滑因子 (FFT平滑实验数据): {smooth_factor}")
    print(f"  Experimental data file: {exp_file}")
    print(f"  实验数据文件: {exp_file}")
    
    print(f"\nUsing parameters: sigma={sigma}, shift={shift}, scale_factor={scale_factor}, smooth_factor={smooth_factor}")
    
    # Check if folders exist
    # 检查文件夹是否存在
    if not os.path.exists(opt_dir):
        print(f"Error: Folder does not exist: {opt_dir}")
        print(f"错误: 文件夹不存在: {opt_dir}")
        print(f"Current directory: {os.getcwd()}")
        print(f"当前目录: {os.getcwd()}")
        print(f"Directory contents: {os.listdir('.')}")
        print(f"目录内容: {os.listdir('.')}")
        return
    
    if not os.path.exists(ecd_dir):
        print(f"Error: Folder does not exist: {ecd_dir}")
        print(f"错误: 文件夹不存在: {ecd_dir}")
        print(f"Directory contents: {os.listdir('.')}")
        print(f"目录内容: {os.listdir('.')}")
        return
    
    # Calculate spectrum
    # 计算谱图
    energy_grid, spectrum = calculate_weighted_spectrum(opt_dir, ecd_dir, sigma)
    
    if energy_grid is None:
        print("Calculation failed")
        print("计算失败")
        return
    
    # ===== Spectrum shift function (units: eV) =====
    # ===== 光谱平移功能（单位：eV）=====
    energy_grid = energy_grid + shift
    # ==================================
    
    # Only keep data in 200-400nm range
    # 只保留200-400nm范围内的数据
    lambda_min = 200
    lambda_max = 400
    
    wavelength_grid = 1239.84 / energy_grid
    mask = (wavelength_grid >= lambda_min) & (wavelength_grid <= lambda_max)
    energy_grid = energy_grid[mask]
    spectrum = spectrum[mask]
    wavelength_grid = wavelength_grid[mask]
    
    # Read experimental data
    # 读取实验数据
    exp_wavelengths, exp_intensities = read_experimental_data(exp_file)

    # FFT smoothing 
    # FFT 平滑 
    if exp_wavelengths is not None and exp_intensities is not None:
        exp_intensities = fft_smooth(exp_intensities, smooth_factor)
        print("Experimental data smoothing completed")
        print("实验数据平滑完成")
    else:
        print("Warning: No experimental data read, will only process calculated data")
        print("警告: 未读取到实验数据，将只处理计算数据")
    
    # Save combined CSV file
    # 保存合并的CSV文件
    if exp_wavelengths is not None and exp_intensities is not None:
        save_combined_spectrum_csv(wavelength_grid, spectrum, exp_wavelengths, exp_intensities, sigma, shift, scale_factor, smooth_factor)
    else:
        print("Warning: No experimental data read, only saving calculated data")
        print("警告: 未读取到实验数据，只保存计算数据")
        # Even without experimental data, save calculated data
        # 即使没有实验数据，也保存计算数据
        filename = f'ecd_calc_only_sigma{sigma}_shift{shift}_scale{scale_factor}.csv'
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("# Calculated ECD spectrum data\n")
            f.write("# 计算ECD谱数据\n")
            f.write("# sigma = {} eV, shift = {} eV, scale_factor = {}\n".format(sigma, shift, scale_factor))
            f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted\n")
            
            # Normalize and scale calculated data
            # 归一化并缩放计算数据
            max_calc_intensity = np.max(np.abs(spectrum))
            if max_calc_intensity == 0:
                norm_calc_spectrum = spectrum
            else:
                norm_calc_spectrum = spectrum / max_calc_intensity
            scaled_calc_spectrum = norm_calc_spectrum * scale_factor
            inverted_calc_spectrum = -scaled_calc_spectrum
            
            for wl, scaled_inten, inverted_inten in zip(wavelength_grid, scaled_calc_spectrum, inverted_calc_spectrum):
                f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f}\n")
        print(f"Calculated data saved to {filename}")
        print(f"计算数据已保存到 {filename}")
    
    # Plot spectrum (with experimental data comparison if available)
    # 绘制谱图（包含实验数据对比）
    plot_spectrum(energy_grid, spectrum, sigma, shift, scale_factor, exp_wavelengths, exp_intensities, plot_exp_name, plot_calc_name, plot_calc_name_invert, smooth_factor)
    
    print("="*60)
    print("Processing completed!")
    print("处理完成!")
    print("="*60)
    
    # Display parameter suggestions
    # 显示参数建议
    if exp_wavelengths is not None and exp_intensities is not None:
        print("\nTip: If calculated data doesn't match experimental data well, adjust scale_factor parameter")
        print("提示: 如果计算数据与实验数据高度不匹配，可以调整scale_factor参数")
        print("For example, set scale_factor to the ratio of experimental data max intensity to calculated data max intensity")
        print("例如，将scale_factor设置为实验数据最大强度与计算数据最大强度的比值")

if __name__ == "__main__":
    main()
