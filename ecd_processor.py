#!/usr/bin/env python3
"""
ORCA ECD数据处理脚本 - 最终版
增加缩放因子和倒置数据功能
"""

import os
import re
import numpy as np

# 尝试导入matplotlib，如果失败则只生成数据文件
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("警告: matplotlib未安装，将只生成数据文件")

def read_weights(weights_file):
    """读取权重文件，尝试多种编码"""
    weights = {}
    
    if not os.path.exists(weights_file):
        print(f"错误: 权重文件不存在: {weights_file}")
        return weights
    
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
            print(f"使用 {encoding} 编码成功读取权重文件")
            break
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(f"使用 {encoding} 编码读取时出错: {e}")
            continue
    
    return weights

def parse_ecd_file(filepath, debug=False):
    """鲁棒解析 ORCA ECD 输出（兼容多种列格式），返回 (energies_array, R_array) 或 None"""
    import math

    if not os.path.exists(filepath):
        if debug: print(f"parse_ecd_file: 文件不存在: {filepath}")
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
        print(f"警告: 无法以常见编码读取 {filepath}")
        return None

    # 查找 ECD/CD 块的文本范围
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
            print(f"parse_ecd_file: 在 {os.path.basename(filepath)} 中找不到 CD SPECTRUM 标记 (尝试的编码: {used_encoding})")
            print("=== DEBUG START (file head 2000 chars) ===")
            print(content[:2000])
            print("=== DEBUG END ===")
        return None

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

    # 按行筛出包含"->"的过渡行（一般为实际数据行）
    lines = block.splitlines()
    data_lines = []
    for ln in lines:
        if '->' in ln and re.search(r'\d', ln):
            # 排除表头线（含 'Energy', 'Wavelength' 等）和短行
            if len(ln.strip()) < 10:
                continue
            # ignore lines that are clearly separators
            if re.match(r'^\s*-+\s*$', ln):
                continue
            data_lines.append(ln)

    if not data_lines:
        if debug:
            print(f"parse_ecd_file: 未在块中发现数据行 (文件 {os.path.basename(filepath)})")
            print("=== DEBUG BLOCK START ===")
            print(block[:2000])
            print("=== DEBUG BLOCK END ===")
        return None

    energies = []
    R_values = []

    for ln in data_lines:
        # 提取所有可以转换为 float 的 token
        toks = ln.split()
        floats = []
        for t in toks:
            try:
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

        # 如果没有检测到 floats，则尝试从 line 中提取带小数的子串
        if not floats:
            found = re.findall(r'-?\d+\.\d+(?:[eE][\+\-]?\d+)?', ln)
            floats = [float(x) for x in found]

        if len(floats) < 2:
            # 无效行，跳过
            if debug:
                print(f"parse_ecd_file: 跳过行（float不足）: {ln}")
            continue

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

        # R 的位置通常是 the 4th float (index 3) in ORCA outputs above; 若不存在则取第一个负/较大值
        R = None
        if len(floats) >= 4:
            R = floats[3]
        else:
            # 回退策略：找第一个绝对值较大的 float（排除能量/波数/波长）
            cand = [f for f in floats if abs(f) > 1e-3 and (f < -0.5 or f > 0.5)]
            if cand:
                R = cand[0]
            else:
                # 最后退化取 floats[-(number)] 末尾靠近的项
                R = floats[-4] if len(floats) >= 4 else floats[-1]

        # 验证 energy 与 R 看起来合理（energy 正，R 可正可负但不为 NaN）
        try:
            energy = float(energy)
            R = float(R)
            energies.append(energy)
            R_values.append(R)
        except Exception:
            if debug:
                print(f"parse_ecd_file: 转换失败 ln: {ln} floats: {floats}")

    if not energies:
        if debug:
            print(f"parse_ecd_file: 解析后无有效能量 (file {os.path.basename(filepath)})")
        return None

    if debug:
        print(f"parse_ecd_file: 在 {os.path.basename(filepath)} 中解析到 {len(energies)} 个跃迁 (encoding={used_encoding})")

    return np.array(energies), np.array(R_values)

def calculate_weighted_spectrum(opt_dir, ecd_dir, sigma=0.3):
    """计算加权ECD谱"""
    print("开始计算加权ECD谱...")
    
    # 尝试读取权重文件
    weights_file = os.path.join(opt_dir, 'ecd_weights.txt')
    weights = read_weights(weights_file)
    
    if not weights:
        print("错误: 无法获取权重数据")
        return None, None
    
    print(f"读取到 {len(weights)} 个构象的权重")
    
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
            print(f"警告: 无法解析构象 {conf_id} 的数据")
    
    if not all_data:
        print("错误: 没有可用的ECD数据")
        return None, None
    
    print(f"成功解析 {len(all_data)} 个构象的ECD数据")
    
    # 创建能量网格 (1.5-8.5 eV，对应约146-827 nm)
    energy_grid = np.linspace(1.5, 8.5, 1000)
    spectrum = np.zeros_like(energy_grid)
    
    # 计算加权平均谱
    total_weight = sum(item['weight'] for item in all_data.values())
    
    for conf_id, conf_data in all_data.items():
        weight = conf_data['weight'] / total_weight
        
        for energy, R in zip(conf_data['energy'], conf_data['R']):
            # 高斯展宽
            gaussian = np.exp(-(energy_grid - energy)**2 / (2 * sigma**2))
            gaussian /= np.sqrt(2 * np.pi * sigma**2)
            spectrum += weight * R * gaussian
    
    return energy_grid, spectrum

def fft_smooth(y, smooth_factor):
    """
    使用 FFT 低通滤波平滑数据
    smooth_factor: 0~1，越小越平滑，推荐 0.02~0.2
    """
    y = np.array(y)
    n = len(y)

    # FFT
    Y = np.fft.rfft(y)

    # 频率裁剪（低通）
    cutoff = int(len(Y) * smooth_factor)   # 保留前百分之 smooth_factor 的频率
    Y[cutoff:] = 0                         # 高频置零

    # 逆 FFT
    y_smooth = np.fft.irfft(Y, n=n)
    return y_smooth

def read_experimental_data(exp_file='result.csv'):
    """读取实验数据文件（支持两种格式：带XYDATA标记的格式和简单CSV格式）"""
    if not os.path.exists(exp_file):
        print(f"警告: 实验数据文件不存在: {exp_file}")
        return None, None
    
    wavelengths = []
    intensities = []
    
    # 尝试多种编码
    encodings = ['utf-8', 'gbk', 'gb2312', 'latin-1']
    
    for encoding in encodings:
        try:
            with open(exp_file, 'r', encoding=encoding) as f:
                lines = f.readlines()
            
            # 检查文件是否包含XYDATA标记
            has_xydata = any('XYDATA' in line for line in lines)
            
            if has_xydata:
                # 格式1: 包含XYDATA标记的格式
                # 寻找数据开始的位置
                start_idx = -1
                for i, line in enumerate(lines):
                    if 'XYDATA' in line:
                        start_idx = i + 1
                        break
                
                if start_idx == -1:
                    print("警告: 实验数据文件中未找到XYDATA标记")
                    return None, None
                
                # 读取数据
                for line in lines[start_idx:]:
                    # 跳过空行和注释行
                    if not line.strip() or ',' not in line:
                        continue
                    
                    # 检查是否到达注释部分
                    if '##### Extended Information' in line or '[Comments]' in line:
                        break
                    
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
                # 格式2: 简单CSV格式（每行: 波长,强度）
                for line in lines:
                    # 跳过空行和注释行（以#开头的行）
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # 尝试解析数据行
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            # 处理科学计数法（如2.18646E-05）
                            wavelength_str = parts[0].strip()
                            intensity_str = parts[1].strip()
                            
                            # 替换科学计数法中的E/e为python可识别的格式
                            if 'E' in intensity_str or 'e' in intensity_str:
                                # 直接使用float()，它能处理科学计数法
                                intensity = float(intensity_str)
                            else:
                                intensity = float(intensity_str)
                            
                            wavelength = float(wavelength_str)
                            
                            wavelengths.append(wavelength)
                            intensities.append(intensity)
                        except ValueError as e:
                            print(f"解析行时出错 '{line}': {e}")
                            continue
            
            if wavelengths:  # 如果读取到数据
                print(f"使用 {encoding} 编码成功读取实验数据文件")
                break
                
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(f"使用 {encoding} 编码读取实验数据时出错: {e}")
            continue
    
    if not wavelengths:
        print("警告: 未读取到实验数据")
        return None, None
    
    # 转换为numpy数组并确保顺序正确（从低波长到高波长）
    wavelengths = np.array(wavelengths)
    intensities = np.array(intensities)
    
    # 检查数据顺序：如果波长是降序（如400, 399.9, 399.8...），则反转
    if len(wavelengths) > 1 and wavelengths[0] > wavelengths[-1]:
        print("检测到波长降序排列，正在反转数据...")
        wavelengths = wavelengths[::-1]
        intensities = intensities[::-1]
    
    print(f"读取到 {len(wavelengths)} 个实验数据点")
    print(f"波长范围: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} nm")
    print(f"强度范围: {intensities.min():.6f} - {intensities.max():.6f}")
    
    return wavelengths, intensities

def save_combined_spectrum_csv(calc_wavelength_grid, calc_spectrum, exp_wavelengths, exp_intensities, sigma, shift, scale_factor, smooth_factor=None):
    """保存合并的谱图数据到CSV，包含计算数据（归一化并缩放）和实验数据（平滑后）"""
    filename = f'ecd_combined_spectrum_sigma{sigma}_shift{shift}_scale{scale_factor}_smooth{smooth_factor}.csv'
    print(f"保存合并的CSV数据到 {filename}...")
    
    # ===== 计算数据的归一化 =====
    max_calc_intensity = np.max(np.abs(calc_spectrum))
    if max_calc_intensity == 0:
        norm_calc_spectrum = calc_spectrum
    else:
        norm_calc_spectrum = calc_spectrum / max_calc_intensity
    
    # ===== 应用缩放因子 =====
    scaled_calc_spectrum = norm_calc_spectrum * scale_factor
    
    # ===== 生成倒置数据 =====
    inverted_calc_spectrum = -scaled_calc_spectrum
    # ======================
    
    # 写入CSV文件
    with open(filename, 'w', encoding='utf-8') as f:
        # 写入文件头
        f.write("# ECD谱数据 - 计算数据（归一化并缩放）与实验数据（平滑）对比\n")
        f.write("# 注意：计算数据已归一化（最大强度 = 1）并乘以缩放因子{}\n".format(scale_factor))
        f.write("# 平滑后实验数据（mdeg单位）\n")
        f.write("# sigma = {} eV, shift = {} eV, scale_factor = {}, smooth_factor = {}\n".format(sigma, shift, scale_factor, smooth_factor))
        f.write("\n")
        
        # 写入计算数据标题
        f.write("=== 计算数据（归一化并缩放，及倒置版本）===\n")
        f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted\n")
        
        # 写入计算数据
        for wl, scaled_inten, inverted_inten in zip(calc_wavelength_grid, scaled_calc_spectrum, inverted_calc_spectrum):
            f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f}\n")
        
        f.write("\n")
        
        # 写入实验数据标题
        f.write("=== 实验数据（平滑）===\n")
        f.write("Wavelength(nm),Experimental_Intensity(mdeg)\n")
        
        # 写入实验数据
        for wl, inten in zip(exp_wavelengths, exp_intensities):
            f.write(f"{wl:.4f},{inten:.6f}\n")
    
    print(f"合并的CSV数据保存完成: {filename}")
    
    # 同时保存一个Origin友好格式（对齐的列）
    origin_filename = f'ecd_origin_format_sigma{sigma}_shift{shift}_scale{scale_factor}_smooth{smooth_factor}.csv'
    print(f"保存Origin格式CSV数据到 {origin_filename}...")
    
    # 创建实验数据的插值函数（线性插值）
    try:
        # 检查是否有足够的实验数据点进行插值
        if len(exp_wavelengths) > 1:
            # 手动实现线性插值（避免依赖scipy）
            exp_interp_intensities = np.interp(calc_wavelength_grid, exp_wavelengths, exp_intensities, 
                                               left=np.nan, right=np.nan)
            
            # 写入Origin友好格式
            with open(origin_filename, 'w', encoding='utf-8') as f:
                f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted,Experimental_Intensity(mdeg)\n")
                for i, wl in enumerate(calc_wavelength_grid):
                    scaled_inten = scaled_calc_spectrum[i]
                    inverted_inten = inverted_calc_spectrum[i]
                    exp_inten = exp_interp_intensities[i]
                    f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f},{exp_inten:.6f}\n")
            
            print(f"Origin格式CSV数据保存完成: {origin_filename}")
        else:
            print("实验数据点太少，无法进行插值，跳过Origin格式保存")
    except Exception as e:
        print(f"插值失败，跳过Origin格式保存: {e}")
        
        # 保存简化版本
        with open(origin_filename, 'w', encoding='utf-8') as f:
            f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted\n")
            for i, wl in enumerate(calc_wavelength_grid):
                scaled_inten = scaled_calc_spectrum[i]
                inverted_inten = inverted_calc_spectrum[i]
                f.write(f"{wl:.4f},{scaled_inten:.6f},{inverted_inten:.6f}\n")
        print("已保存仅含计算数据的简化版本")

def plot_spectrum(energy_grid, calc_spectrum, sigma, shift, scale_factor, exp_wavelengths=None, exp_intensities=None, plot_exp_name=None, plot_calc_name=None, plot_calc_name_invert=None, smooth_factor=None):
    if not HAS_MATPLOTLIB:
        print("matplotlib未安装，跳过绘图")
        return
    
    # 设置Times New Roman字体
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'  # 数学字体也设置为Times New Roman风格

    # ===== 计算数据的归一化和缩放 =====
    max_calc_intensity = np.max(np.abs(calc_spectrum))
    if max_calc_intensity == 0:
        norm_calc_spectrum = calc_spectrum
    else:
        norm_calc_spectrum = calc_spectrum / max_calc_intensity
    
    # 应用缩放因子
    scaled_calc_spectrum = norm_calc_spectrum * scale_factor
    
    # 生成倒置数据
    inverted_calc_spectrum = -scaled_calc_spectrum
    # ======================
    
    plt.figure(figsize=(10, 6))
    
    # 计算数据：主坐标改为波长
    wavelength_grid = 1239.84 / energy_grid
    
    # 绘制计算谱图（缩放后）
    plt.plot(wavelength_grid, scaled_calc_spectrum, 'r--', linewidth=1.5, 
             label=f'{plot_calc_name}')
    
    # 绘制倒置计算谱图
    plt.plot(wavelength_grid, inverted_calc_spectrum, 'b--', linewidth=1.5, alpha=0.7,
             label=f'{plot_calc_name_invert}')
    
    # 绘制实验数据（如果提供）
    if exp_wavelengths is not None and exp_intensities is not None:
        # 实验数据不进行归一化，保持原始值
        plt.plot(exp_wavelengths, exp_intensities, 'k-', linewidth=1.5, label=f'{plot_exp_name}')
    
    # 主 x 轴：200–400 nm，均匀刻度
    plt.xlim(200, 400)
    plt.xticks([200, 250, 300, 350, 400])
    plt.xlabel('Wavelength (nm)', fontsize=14, fontweight='bold')
    
    # y 轴
    plt.ylabel('Δε (M$^{-1}$cm$^{-1}$)', fontsize=14, fontweight='bold')
    
    # 根据是否有实验数据设置标题
    if exp_wavelengths is not None and exp_intensities is not None:
        plt.title(f'ECD Spectrum Comparison\nσ={sigma} eV, shift={shift} eV, scale={scale_factor}, smooth={smooth_factor}', fontsize=14)
    else:
        plt.title(f'Calculated ECD Spectrum\nσ={sigma} eV, shift={shift} eV, scale={scale_factor}', fontsize=14)
    
    # 零线
    plt.axhline(y=0, color='black', linewidth=0.5)
    
    # 副 x 轴（上方）：能量 eV
    ax = plt.gca()
    
    def nm_to_ev(x):
        return 1239.84 / x
    
    def ev_to_nm(x):
        return 1239.84 / x
    
    ax2 = ax.secondary_xaxis('top', functions=(nm_to_ev, ev_to_nm))
    ax2.set_xlabel('Energy (eV)', fontsize=12)
    
    # 添加图例
    plt.legend(loc='best', shadow=True, edgecolor='black', facecolor='white', fontsize=10)
    
    # 网格
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # 根据是否有实验数据设置输出文件名
    if exp_wavelengths is not None and exp_intensities is not None:
        output_file = f'ecd_spectrum_comparison_sigma{sigma}_shift{shift}_scale{scale_factor}_smooth{smooth_factor}.png'
    else:
        output_file = f'ecd_spectrum_sigma{sigma}_shift{shift}_scale{scale_factor}.png'
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"谱图已保存为 {output_file}")
    plt.show()

def main():
    """主函数"""
    print("="*60)
    print("ORCA ECD数据处理脚本 - 增强版")
    print("="*60)
    
    # 设置计算参数
    opt_dir = 'opt_conf'    # 优化结果文件夹
    ecd_dir = 'ecd_opt'     # ECD结果文件夹
    sigma = 0.2             # 高斯展宽参数
    shift = 0.18             # 能谱平移（负值向右移，如-0.2）
    scale_factor = 8        # 计算数据缩放因子（用于与实验数据对齐）

    # 设置实验参数
    exp_file = '74-1-50-1110-P.csv' # 实验数据文件
    smooth_factor = 0.03    #平滑实验数据的FFT因子（0~1，越小越平滑）

    # 设置图例名称
    plot_exp_name = 'Experimental ECD of A13B5'
    plot_calc_name = "Calculated ECD of (7$R$, 8$R$, 7'$S$, 8'$R$)- A13B5"
    plot_calc_name_invert = "Calculated ECD of (7$S$, 8$S$, 7'$R$, 8'$S$)- A13B5"
    
    # 从用户获取参数（可选）
    print("当前参数设置：")
    print(f"  sigma (高斯展宽): {sigma}")
    print(f"  shift (能谱平移): {shift}")
    print(f"  scale_factor (计算数据缩放因子): {scale_factor}")
    print(f"  平滑因子 (FFT平滑实验数据): {smooth_factor}")
    print(f"  实验数据文件: {exp_file}")
    
    print(f"\n使用参数: sigma={sigma}, shift={shift}, scale_factor={scale_factor}, smooth_factor={smooth_factor}")
    
    # 检查文件夹是否存在
    if not os.path.exists(opt_dir):
        print(f"错误: 文件夹不存在: {opt_dir}")
        print(f"当前目录: {os.getcwd()}")
        print(f"目录内容: {os.listdir('.')}")
        return
    
    if not os.path.exists(ecd_dir):
        print(f"错误: 文件夹不存在: {ecd_dir}")
        print(f"目录内容: {os.listdir('.')}")
        return
    
    # 计算谱图
    energy_grid, spectrum = calculate_weighted_spectrum(opt_dir, ecd_dir, sigma)
    
    if energy_grid is None:
        print("计算失败")
        return
    
    # ===== 光谱平移功能（单位：eV）=====
    energy_grid = energy_grid + shift
    # ==================================
    
    # 只保留200-400nm范围内的数据
    lambda_min = 200
    lambda_max = 400
    
    wavelength_grid = 1239.84 / energy_grid
    mask = (wavelength_grid >= lambda_min) & (wavelength_grid <= lambda_max)
    energy_grid = energy_grid[mask]
    spectrum = spectrum[mask]
    wavelength_grid = wavelength_grid[mask]
    
    # 读取实验数据
    exp_wavelengths, exp_intensities = read_experimental_data(exp_file)

    # FFT 平滑 
    if exp_wavelengths is not None and exp_intensities is not None:
        exp_intensities = fft_smooth(exp_intensities, smooth_factor)
        print("实验数据平滑完成")
    else:
        print("警告: 未读取到实验数据，将只处理计算数据")
    
    # 保存合并的CSV文件
    if exp_wavelengths is not None and exp_intensities is not None:
        save_combined_spectrum_csv(wavelength_grid, spectrum, exp_wavelengths, exp_intensities, sigma, shift, scale_factor, smooth_factor)
    else:
        print("警告: 未读取到实验数据，只保存计算数据")
        # 即使没有实验数据，也保存计算数据
        filename = f'ecd_calc_only_sigma{sigma}_shift{shift}_scale{scale_factor}.csv'
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("# 计算ECD谱数据\n")
            f.write("# sigma = {} eV, shift = {} eV, scale_factor = {}\n".format(sigma, shift, scale_factor))
            f.write("Wavelength(nm),Calculated_Scaled,Calculated_Inverted\n")
            
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
        print(f"计算数据已保存到 {filename}")
    
    # 绘制谱图（包含实验数据对比）
    plot_spectrum(energy_grid, spectrum, sigma, shift, scale_factor, exp_wavelengths, exp_intensities, plot_exp_name, plot_calc_name, plot_calc_name_invert, smooth_factor)
    
    print("="*60)
    print("处理完成!")
    print("="*60)
    
    # 显示参数建议
    if exp_wavelengths is not None and exp_intensities is not None:
        print("\n提示: 如果计算数据与实验数据高度不匹配，可以调整scale_factor参数")
        print("例如，将scale_factor设置为实验数据最大强度与计算数据最大强度的比值")

if __name__ == "__main__":
    main()