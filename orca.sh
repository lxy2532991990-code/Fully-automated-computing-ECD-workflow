#!/bin/bash
#SBATCH -J orca-70-ecd
#SBATCH -p cn-short
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o orca-70-calc_%j.out
#SBATCH -e orca-70-calc_%j.err
#SBATCH -A jiangy_g1
#SBATCH --qos=jiangycns
#SBATCH --no-requeue

# 清除环境
module purge

# 设置 ORCA 环境变量
export ORCA_DIR=/home/jiangy/jiangy_liuxy/lustre1/orca-6.1.0-f.0_linux_x86-64
export PATH=$ORCA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$ORCA_DIR/lib:$LD_LIBRARY_PATH

# 设置 OpenMPI
export MPI_HOME=/home/jiangy/jiangy_liuxy/lustre1/openmpi-4.1.8-install
export PATH=$MPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1

cd /home/jiangy/jiangy_liuxy/lustre1/ECD/B9D3A-RSSR-TEST

echo "=== 环境验证 ==="
echo "ORCA路径: $ORCA_DIR"
echo "MPI路径: $MPI_HOME"
echo "PATH: $PATH"
which orca
which mpirun

echo "=== 开始构象搜索 ==="

# 构象搜索输入文件
conf_search_base="conf-70"  # 构象搜索的基本文件名
conf_search_inp="${conf_search_base}.inp"
initial_xyz="70.xyz"  # 初始构象文件

# 检查初始构象文件是否存在
if [ ! -f "$initial_xyz" ]; then
    echo "错误：初始构象文件 $initial_xyz 不存在"
    exit 1
fi

# 创建构象搜索输入文件
cat > "$conf_search_inp" << EOF
!GOAT XTB

%pal nprocs 64 end
%maxcore 6000

* xyzfile 0 1 $initial_xyz

EOF

echo "=== 执行构象搜索 ==="
conf_search_dir="conformational_search"
mkdir -p $conf_search_dir
$ORCA_DIR/bin/orca "$conf_search_inp" > "${conf_search_dir}/${conf_search_base}.out"

# 检查构象搜索是否成功
if [ $? -ne 0 ]; then
    echo "错误：构象搜索失败"
    exit 1
fi

echo "构象搜索完成"

# 查找构象搜索输出文件（根据ORCA命名规则）
conformer_output="${conf_search_base}.finalensemble.xyz"

# 检查文件是否存在
if [ ! -f "$conformer_output" ]; then
    echo "错误：构象搜索结果文件 $conformer_output 不存在"
    echo "正在查找可能的输出文件..."
    ls -la *.finalensemble* *.xyz *.out 2>/dev/null
    exit 1
fi

echo "找到构象搜索结果文件: $conformer_output"

# 复制文件为多构象文件
multi_xyz_file="multi_conformers.xyz"
cp "$conformer_output" "$multi_xyz_file"

echo "=== 开始拆分构象文件 ==="
# 设置参数
atoms_per_conformer=59  # 每个构象的原子数
max_conformers_to_process=30  # 最多处理的构象数量

# 拆分多构象文件为单个构象文件
conformer_count=0
temp_file=$(mktemp)

# 使用awk来拆分文件
awk -v atoms="$atoms_per_conformer" '
BEGIN { count = 0 }
{
    if (NR % (atoms + 2) == 1) {
        # 新构象开始
        if (count > 0) close("temp_conformer_" count-1 ".xyz")
        count++
        file = "temp_conformer_" count-1 ".xyz"
    }
    print > file
}
END { print count }
' "$multi_xyz_file" > "$temp_file"

conformer_count=$(cat "$temp_file")
rm "$temp_file"

echo "检测到 $conformer_count 个构象"

# 确定实际要处理的构象数量
if [ $conformer_count -gt $max_conformers_to_process ]; then
    echo "只处理前 $max_conformers_to_process 个构象（总共 $conformer_count 个）"
    conformers_to_process=$max_conformers_to_process
else
    conformers_to_process=$conformer_count
    echo "处理所有 $conformer_count 个构象"
fi

# 重命名并保留需要的构象文件
for ((i=0; i<conformers_to_process; i++)); do
    if [ -f "temp_conformer_${i}.xyz" ]; then
        mv "temp_conformer_${i}.xyz" "conf-${i}.xyz"
    fi
done

# 清理不需要的临时构象文件
for ((i=conformers_to_process; i<conformer_count; i++)); do
    rm -f "temp_conformer_${i}.xyz" 2>/dev/null
done

# 设置临时目录到有足够空间的地方
export TMPDIR="/home/jiangy/jiangy_liuxy/lustre1/temp"
mkdir -p $TMPDIR

echo "=== 创建单个构象计算脚本 ==="

# 创建单个构象计算脚本
cat > "single_conformer_calc.sh" << 'EOF'
#!/bin/bash
#SBATCH -J orca-conf
#SBATCH -p cn-short
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o orca-conf_%j.out
#SBATCH -e orca-conf_%j.err
#SBATCH -A jiangy_g1
#SBATCH --qos=jiangycns
#SBATCH --no-requeue

log_suffix="conf-${1}"
exec > >(tee "orca-${log_suffix}.out") 2> >(tee "orca-${log_suffix}.err" >&2)

# 清除环境
module purge

# 设置 ORCA 环境变量
export ORCA_DIR=/home/jiangy/jiangy_liuxy/lustre1/orca-6.1.0-f.0_linux_x86-64
export PATH=$ORCA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$ORCA_DIR/lib:$LD_LIBRARY_PATH

# 设置 OpenMPI
export MPI_HOME=/home/jiangy/jiangy_liuxy/lustre1/openmpi-4.1.8-install
export PATH=$MPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
export TMPDIR="/home/jiangy/jiangy_liuxy/lustre1/temp"

cd /home/jiangy/jiangy_liuxy/lustre1/ECD/B9D3A-RSSR-TEST

i=$1
echo "=== 处理构象 $i ==="
echo "开始时间: $(date)"

# 第一步：几何优化
opt_inp_file="opt-conf-${i}.inp"
xyz_file="conf-${i}.xyz"

# 检查xyz文件是否存在
if [ ! -f "$xyz_file" ]; then
    echo "错误：文件 $xyz_file 不存在"
    exit 1
fi

echo "=== 构象 $i: 几何优化 ==="

# 创建几何优化输入文件内容
cat > "$opt_inp_file" << EOL
! B3LYP OPT FREQ DEF2-TZVP D4 CPCM(Methanol) 
! TightOpt
%pal nprocs 64 end
%maxcore 6000

%cpcm
end

%freq
   Temp 298.15
end

* xyzfile 0 1 $xyz_file

EOL

# 处理几何优化文件
opt_base_name="${opt_inp_file%.inp}"
echo "开始几何优化计算构象 $i..."
$ORCA_DIR/bin/orca "$opt_inp_file" > "${opt_base_name}.out" 2>&1
opt_exit_code=$?

# 检查几何优化计算是否成功
if [ $opt_exit_code -eq 0 ]; then
    echo "构象 $i 几何优化完成"
else
    echo "错误：构象 $i 几何优化失败，退出码: $opt_exit_code"
fi

# 第二步：ECD计算
echo "=== 构象 $i: ECD计算 ==="

# 优化后的结构文件
opt_xyz_file="opt-conf-${i}.xyz"
ECD_inp_file="ECD-conf-${i}.inp"

# 检查优化后的结构文件是否存在
if [ ! -f "$opt_xyz_file" ]; then
    echo "错误：优化后的结构文件 $opt_xyz_file 不存在"
    # 如果几何优化失败，跳过ECD计算
    echo "跳过构象 $i 的ECD计算"
    exit 0
fi

# 创建ECD计算输入文件内容
cat > "$ECD_inp_file" << EOL
!B3LYP DEF2-TZVP CPCM(METHANOL)

%pal nprocs 64 end
%maxcore 6000

%TDDFT
      NROOTS  25
      TDA     FALSE
END

* xyzfile 0 1 opt-conf-${i}.xyz

EOL

# 处理ECD计算文件
ECD_base_name="${ECD_inp_file%.inp}"
echo "开始ECD计算构象 $i..."
$ORCA_DIR/bin/orca "$ECD_inp_file" > "${ECD_base_name}.out" 2>&1
ECD_exit_code=$?

# 检查ECD计算是否成功
if [ $ECD_exit_code -eq 0 ]; then
    echo "构象 $i ECD计算完成"
else
    echo "错误：构象 $i ECD计算失败，退出码: $ECD_exit_code"
fi

# 立即清理临时文件以释放空间
find . -name "${opt_base_name}.*.tmp" -delete 2>/dev/null
find . -name "${ECD_base_name}.*.tmp" -delete 2>/dev/null

echo "=== 构象 $i 处理完成 ==="
echo "结束时间: $(date)"
EOF

chmod +x single_conformer_calc.sh

echo "=== 创建任务池管理函数 ==="

# 函数：检查文件是否存在并记录
check_and_move_files() {
    local i=$1
    local opt_conf_dir=$2
    local ECD_opt_dir=$3
    
    echo "移动构象 $i 的输出文件..."
    
    # 检查并移动几何优化相关文件
    declare -a opt_files=(
        "opt-conf-${i}.inp"
        "opt-conf-${i}.xyz" 
        "opt-conf-${i}.out"
    )
    
    for file in "${opt_files[@]}"; do
        if ls $file 1> /dev/null 2>&1; then
            echo "移动文件: $file -> $opt_conf_dir/"
            mv $file ${opt_conf_dir}/ 2>/dev/null || echo "移动 $file 失败"
        fi
    done
    
    # 检查并移动ECD计算相关文件
    declare -a ECD_files=(
        "ECD-conf-${i}.inp"
        "ECD-conf-${i}.out"
    )
    
    for file in "${ECD_files[@]}"; do
        if ls $file 1> /dev/null 2>&1; then
            echo "移动文件: $file -> $ECD_opt_dir/"
            mv $file ${ECD_opt_dir}/ 2>/dev/null || echo "移动 $file 失败"
        fi
    done
    
    # 清理临时文件
    find . -name "opt-conf-${i}.*.tmp" -delete 2>/dev/null
    find . -name "ECD-conf-${i}.*.tmp" -delete 2>/dev/null
}

# 函数：维护任务池，当有任务完成时提交新任务
manage_job_pool() {
    local total_jobs=$1
    local pool_size=$2
    local submitted=0
    local completed=0
    local max_wait_time=86400  # 最大等待时间24小时
    local start_time=$(date +%s)
    
    # 存储运行中的作业ID和对应的构象编号
    declare -A running_jobs
    declare -a job_queue=()
    
    echo "开始任务池管理: 总共 $total_jobs 个任务，池大小 $pool_size"
    
    # 函数：提交单个任务
    submit_single_job() {
        local i=$1
        echo "提交构象 $i 的计算任务..."
        job_output=$(sbatch --parsable single_conformer_calc.sh $i)
        job_id=$job_output
        running_jobs[$job_id]=$i
        job_queue+=($job_id)
        echo "构象 $i 作业ID: $job_id" | tee -a $job_log
        echo "作业 $job_id (构象 $i) 已提交" >> $job_log
    }
    
    # 函数：检查作业状态并更新运行中任务列表
    check_job_status() {
        local new_job_queue=()
        local current_completed=$completed
        
        for job_id in "${job_queue[@]}"; do
            local conf_index=${running_jobs[$job_id]}
            
            if squeue -j $job_id &> /dev/null; then
                # 作业还在运行
                new_job_queue+=($job_id)
            else
                # 作业已完成，检查退出状态
                job_state=$(sacct -j $job_id --format=State -n -P -S $(date -d "2 days ago" +%Y-%m-%d) 2>/dev/null | head -1)
                if [ "$job_state" == "COMPLETED" ]; then
                    echo "作业 $job_id (构象 $conf_index) 已完成"
                    # 移动该构象的输出文件
                    check_and_move_files $conf_index $opt_conf_dir $ECD_opt_dir
                    ((completed++))
                else
                    echo "作业 $job_id (构象 $conf_index) 失败，状态: $job_state"
                    # 即使失败也移动已有文件并计数
                    check_and_move_files $conf_index $opt_conf_dir $ECD_opt_dir
                    ((completed++))
                fi
                # 从运行列表中移除
                unset running_jobs[$job_id]
            fi
        done
        
        # 更新作业队列
        job_queue=("${new_job_queue[@]}")
        
        # 如果有任务完成，打印状态
        if [ $completed -gt $current_completed ]; then
            echo "有 $(($completed - $current_completed)) 个任务完成，当前运行: ${#running_jobs[@]}, 总计完成: $completed/$total_jobs"
        fi
    }
    
    # 初始提交一批任务
    echo "=== 初始提交 $pool_size 个任务 ==="
    for ((i=0; i<pool_size && submitted<total_jobs; i++)); do
        submit_single_job $submitted
        ((submitted++))
        sleep 10  # 短暂间隔避免系统负载
    done
    
    # 主循环：监控任务状态并补充新任务
    while [ $completed -lt $total_jobs ]; do
        # 检查当前时间是否超过最大等待时间
        local current_time=$(date +%s)
        local elapsed_time=$((current_time - start_time))
        
        if [ $elapsed_time -gt $max_wait_time ]; then
            echo "警告：总等待时间超过24小时"
            break
        fi
        
        # 检查作业状态
        check_job_status
        
        # 如果有任务完成且还有未提交的任务，提交新任务
        while [ ${#running_jobs[@]} -lt $pool_size ] && [ $submitted -lt $total_jobs ]; do
            submit_single_job $submitted
            ((submitted++))
            sleep 5
        done
        
        # 如果还有任务在运行，等待一段时间再检查
        if [ ${#running_jobs[@]} -gt 0 ]; then
            sleep 60  # 等待1分钟再检查
        else
            # 没有运行中的任务，检查是否全部完成
            if [ $completed -eq $total_jobs ]; then
                break
            else
                sleep 30
            fi
        fi
    done
	
    echo "=== Final re-scan for completed jobs ==="
    for ((i=0; i<total_jobs; i++)); do
        if [ -f "opt-conf-${i}.out" ] || [ -f "ECD-conf-${i}.out" ]; then
            echo "检测到构象 $i 的输出文件，执行最终移动..."
            check_and_move_files $i $opt_conf_dir $ECD_opt_dir
        fi
    done
    
    echo "所有任务已完成: $completed/$total_jobs"
}

echo "=== 使用任务池方式提交计算任务 ==="

# 任务池参数
pool_size=30  # 同时运行的任务数量
total_jobs=$conformers_to_process  # 总任务数

echo "总共 $total_jobs 个构象，使用任务池管理，同时运行 $pool_size 个任务"

# 创建输出目录
opt_conf_dir="opt_conf"
ECD_opt_dir="ECD_opt"
mkdir -p ${opt_conf_dir}
mkdir -p ${ECD_opt_dir}

# 创建作业日志文件
job_log="job_submission.log"
echo "任务池作业提交日志" > $job_log
echo "==================" >> $job_log
echo "开始时间: $(date)" >> $job_log
echo "总任务数: $total_jobs, 池大小: $pool_size" >> $job_log

# 使用任务池管理函数提交任务
manage_job_pool $total_jobs $pool_size

echo "=== 所有计算任务已完成 ==="
echo "共处理了 $total_jobs 个构象"
echo "几何优化结果保存在: $opt_conf_dir"
echo "ECD计算结果保存在: $ECD_opt_dir"
echo "主作业结束时间: $(date)"
echo "结束时间: $(date)" >> $job_log

# 创建玻尔兹曼分布计算脚本
echo "=== 创建玻尔兹曼分布计算脚本 ==="

cat > "boltzmann_calculator.py" << 'PYEOF'
#!/usr/bin/env python3
import re
import math
import os
import glob
import sys

def parse_gibbs_energy_from_out(filename):
    """从ORCA .out文件中提取吉布斯自由能"""
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        # 如果UTF-8失败，尝试其他编码
        with open(filename, 'r', encoding='latin-1') as f:
            content = f.read()
    
    # 匹配吉布斯自由能行 - 更灵活的正则表达式
    patterns = [
        r"Final Gibbs free energy\s*\.+\s*([-\d.]+)\s*Eh",
        r"G-E\(el\)\s*\.+\s*([-\d.]+)\s*Eh",  # 备选模式
        r"Total Gibbs free energy\s*\.+\s*([-\d.]+)\s*Eh"
    ]
    
    for pattern in patterns:
        match = re.search(pattern, content)
        if match:
            return float(match.group(1))
    
    # 如果上面的模式没找到，尝试更宽松的搜索
    gibbs_pattern = r"Gibbs.*free.*energy.*?([-\d.]+)\s*Eh"
    match = re.search(gibbs_pattern, content, re.IGNORECASE)
    if match:
        return float(match.group(1))
    
    raise ValueError(f"在文件 {filename} 中未找到吉布斯自由能")

def calculate_boltzmann_distribution(energy_list, temperature=298.15):
    """计算玻尔兹曼分布"""
    # 转换常数：Hartree to kcal/mol
    hartree_to_kcal = 627.509
    
    # 转换为kcal/mol
    energies_kcal = [energy * hartree_to_kcal for energy in energy_list]
    
    # 找到最低能量（作为参考）
    min_energy = min(energies_kcal)
    
    # 计算相对能量
    relative_energies = [energy - min_energy for energy in energies_kcal]
    
    # 气体常数 (kcal/mol·K)
    R = 0.001987
    
    # 计算玻尔兹曼因子
    boltzmann_factors = [math.exp(-deltaG / (R * temperature)) 
                        for deltaG in relative_energies]
    
    # 计算配分函数
    partition_function = sum(boltzmann_factors)
    
    # 计算各构象比例
    percentages = [factor / partition_function * 100 for factor in boltzmann_factors]
    
    return percentages, relative_energies, energies_kcal

def main():
    # 使用当前目录下的opt_conf目录
    base_path = "opt_conf"
    
    # 如果通过命令行参数指定了路径，使用该路径
    if len(sys.argv) > 1:
        base_path = sys.argv[1]
    
    # 获取所有.out文件
    out_files = glob.glob(os.path.join(base_path, "*.out"))
    
    if not out_files:
        print(f"在路径 {base_path} 中未找到.out文件")
        print("请确保计算已经完成并且输出文件在正确的位置")
        return
    
    print(f"找到 {len(out_files)} 个.out文件")
    
    gibbs_energies = []
    file_names = []
    conformer_indices = []
    
    print("\n读取吉布斯自由能数据...")
    for file_path in out_files:
        try:
            energy = parse_gibbs_energy_from_out(file_path)
            gibbs_energies.append(energy)
            file_name = os.path.basename(file_path)
            file_names.append(file_name)
            
            # 从文件名中提取构象编号
            match = re.search(r'opt-conf-(\d+)', file_name)
            if match:
                conformer_indices.append(int(match.group(1)))
            else:
                conformer_indices.append(len(gibbs_energies)-1)
                
            print(f"{file_name}: {energy:.8f} Eh")
        except Exception as e:
            print(f"错误读取 {os.path.basename(file_path)}: {e}")
    
    if gibbs_energies:
        percentages, rel_energies, abs_energies = calculate_boltzmann_distribution(gibbs_energies)
        
        print("\n" + "="*80)
        print("玻尔兹曼分布结果 (T = 298.15 K)")
        print("="*80)
        print(f"{'构象编号':<10} {'构象文件':<25} {'绝对能量(Eh)':<15} {'相对能量(kcal/mol)':<20} {'分布比例(%)':<15}")
        print("-"*80)
        
        # 按分布比例排序（从高到低）
        sorted_data = sorted(zip(conformer_indices, file_names, gibbs_energies, rel_energies, percentages), 
                           key=lambda x: x[4], reverse=True)
        
        for i, (conf_idx, filename, abs_energy, rel_energy, percentage) in enumerate(sorted_data):
            print(f"{conf_idx:<10} {filename:<25} {abs_energy:<15.8f} {rel_energy:<20.2f} {percentage:<15.2f}")
        
        print("-"*80)
        print(f"最低能量构象: {min(gibbs_energies):.8f} Eh (构象 {conformer_indices[gibbs_energies.index(min(gibbs_energies))]})")
        print(f"能量跨度: {max(rel_energies):.2f} kcal/mol")
        print(f"主要构象数量(>5%): {sum(1 for p in percentages if p > 5)}")
        print(f"主要构象数量(>1%): {sum(1 for p in percentages if p > 1)}")
        print(f"主要构象比例总和(>1%): {sum(p for p in percentages if p > 1):.2f}%")
        
        # 保存结果到文件
        output_file = os.path.join(base_path, "boltzmann_results.txt")
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("玻尔兹曼分布计算结果\n")
            f.write("="*60 + "\n")
            f.write(f"温度: 298.15 K\n")
            f.write(f"总构象数: {len(out_files)}\n")
            f.write(f"计算时间: {os.popen('date').read().strip()}\n")
            f.write("="*60 + "\n\n")
            f.write(f"{'构象编号':<10} {'构象文件':<25} {'绝对能量(Eh)':<15} {'相对能量(kcal/mol)':<20} {'分布比例(%)':<15}\n")
            f.write("-"*60 + "\n")
            for conf_idx, filename, abs_energy, rel_energy, percentage in sorted_data:
                f.write(f"{conf_idx:<10} {filename:<25} {abs_energy:<15.8f} {rel_energy:<20.2f} {percentage:<15.2f}\n")
            
            f.write("\n" + "="*60 + "\n")
            f.write("汇总信息:\n")
            f.write(f"最低能量构象: {min(gibbs_energies):.8f} Eh (构象 {conformer_indices[gibbs_energies.index(min(gibbs_energies))]})\n")
            f.write(f"能量跨度: {max(rel_energies):.2f} kcal/mol\n")
            f.write(f"主要构象数量(>5%): {sum(1 for p in percentages if p > 5)}\n")
            f.write(f"主要构象数量(>1%): {sum(1 for p in percentages if p > 1)}\n")
            f.write(f"主要构象比例总和(>1%): {sum(p for p in percentages if p > 1):.2f}%\n")
        
        print(f"\n结果已保存到: {output_file}")
        
        # 创建用于ECD谱加权的文件
        ecd_weight_file = os.path.join(base_path, "ecd_weights.txt")
        with open(ecd_weight_file, 'w', encoding='utf-8') as f:
            f.write("# ECD谱加权文件 - 用于加权平均ECD谱\n")
            f.write("# 构象编号, 权重(0-1)\n")
            for conf_idx, filename, abs_energy, rel_energy, percentage in sorted_data:
                weight = percentage / 100.0
                f.write(f"{conf_idx}, {weight:.6f}\n")
        
        print(f"ECD谱加权文件已保存到: {ecd_weight_file}")
        
        # 创建简化的主要构象列表
        main_conformers = [(conf_idx, percentage) for conf_idx, _, _, _, percentage in sorted_data if percentage > 1.0]
        if main_conformers:
            print("\n主要构象列表(>1%):")
            for conf_idx, percentage in main_conformers:
                print(f"  构象 {conf_idx}: {percentage:.2f}%")

if __name__ == "__main__":
    main()
PYEOF

chmod +x boltzmann_calculator.py

# 运行玻尔兹曼分布计算
echo "=== 设置Python环境 ==="
# 加载特定版本的Python
if [ -f "/appsnew/source/Python-3.8.6.sh" ]; then
    source /appsnew/source/Python-3.8.6.sh
    echo "Python版本: $(python --version 2>&1)"
    echo "Python路径: $(which python)"
else
    echo "警告：未找到Python-3.8.6.sh，尝试使用默认python3"
    # 检查python3是否可用
    if command -v python3 &> /dev/null; then
        echo "使用默认python3: $(python3 --version 2>&1)"
    else
        echo "错误：未找到可用的Python解释器"
        echo "请检查Python环境配置"
        exit 1
    fi
fi

echo "=== 开始玻尔兹曼分布计算 ==="
python boltzmann_calculator.py ${opt_conf_dir}

# 检查玻尔兹曼计算是否成功
if [ $? -eq 0 ]; then
    echo "玻尔兹曼分布计算完成"
else
    echo "警告：玻尔兹曼分布计算可能存在问题"
    echo "尝试使用python3..."
    if command -v python3 &> /dev/null; then
        python3 boltzmann_calculator.py ${opt_conf_dir}
    fi
fi

# 创建详细的结果汇总脚本
cat > "detailed_summary.sh" << 'EOF'
#!/bin/bash
echo "=== 详细计算结果汇总 ==="
echo "几何优化完成文件:"
find opt_conf -name "opt-conf-*.xyz" | sort -V | wc -l
find opt_conf -name "opt-conf-*.xyz" | sort -V

echo -e "\nECD计算完成文件:"
find ECD_opt -name "ECD-conf-*.out" | sort -V | wc -l
find ECD_opt -name "ECD-conf-*.out" | sort -V

echo -e "\n玻尔兹曼分布结果:"
if [ -f "opt_conf/boltzmann_results.txt" ]; then
    echo "玻尔兹曼分布结果文件存在"
    echo "主要构象分布:"
    grep -A 5 "主要构象数量" opt_conf/boltzmann_results.txt
else
    echo "玻尔兹曼分布结果文件不存在"
fi

echo -e "\n缺失的几何优化文件:"
for i in {0..49}; do
    if [ ! -f "opt_conf/opt-conf-${i}.xyz" ]; then
        echo "opt-conf-${i}.xyz"
    fi
done

echo -e "\n缺失的ECD计算文件:"
for i in {0..49}; do
    if [ ! -f "ECD_opt/ECD-conf-${i}.out" ]; then
        echo "ECD-conf-${i}.out"
    fi
done

echo -e "\n失败的作业:"
sacct -u $USER --format=JobID,JobName,State,Elapsed,ExitCode -S $(date -d "yesterday" +%Y-%m-%d) | grep -E "FAILED|CANCELLED|TIMEOUT|NODE_FAIL"
EOF

chmod +x detailed_summary.sh

# 创建一个清理脚本
cat > "cleanup.sh" << 'EOF'
#!/bin/bash
echo "=== 清理临时文件 ==="
echo "清理临时Python环境..."
if command -v deactivate &> /dev/null; then
    deactivate
    echo "Python虚拟环境已退出"
fi

echo "清理临时文件..."
# 清理可能的临时文件
find . -name "*.tmp" -delete 2>/dev/null
find . -name "*.temp" -delete 2>/dev/null
find . -name "temp_conformer_*.xyz" -delete 2>/dev/null

echo "临时文件清理完成"
EOF

chmod +x cleanup.sh

echo "可以使用 './detailed_summary.sh' 来查看详细的结果汇总"
echo "玻尔兹曼分布结果保存在: opt_conf/boltzmann_results.txt"
echo "ECD谱加权文件保存在: opt_conf/ecd_weights.txt"
echo "作业提交日志保存在: $job_log"
echo "可以使用 './cleanup.sh' 来清理临时文件"
