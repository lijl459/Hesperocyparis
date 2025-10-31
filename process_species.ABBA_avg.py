import sys
import os

def process_data(input_file, output_file):
    # 创建一个字典来存储每对物种的z值和gamma的总和及计数
    species_pair_data = {}

    # 存储所有出现过的物种
    all_species = set()

    # 检查输入文件是否存在
    if not os.path.isfile(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在!")
        sys.exit(1)

    # 打开输入文件并读取数据
    with open(input_file, 'r') as infile:
        # 跳过第一行（标题行）
        next(infile)
        
        for line in infile:
            # 按制表符拆分每一行
            row = line.strip().split(',')
            if len(row) != 9:
                continue  # 如果某行数据不完整则跳过

            # 读取数据
            species1,species2,species3,species4,ABBA_count,BABA_count,D_value,Z_value,P_value = row
            
            Z_value = float(Z_value)
            D_value = float(D_value)
            ABBA_count = int(ABBA_count)
            BABA_count = int(BABA_count)
            #if ABBA_count + BABA_count < 100 and species4 != "Cusem":
            #    print(line)
            #    continue

            if Z_value < 0 :
                pair = tuple(sorted([species1, species3]))
            else:
                pair = tuple(sorted([species2, species3]))

            abs_z = abs(Z_value)
            abs_d = abs(D_value)
            if abs_z > 1000:
                print(line)

            # 更新字典中的z值和gamma的总和与计数，用于后续取平均
            if pair in species_pair_data:
                existing_z_sum, existing_d_sum, count = species_pair_data[pair]
                species_pair_data[pair] = (existing_z_sum + abs_z, existing_d_sum + abs_d, count + 1)
            else:
                species_pair_data[pair] = (abs_z, abs_d, 1)





    # 将结果按物种对（sp1, sp2）字典序排序
    sorted_pairs = sorted(species_pair_data.items(), key=lambda item: (item[0][0], item[0][1]))

    # 将结果写入输出文件
    with open(output_file, 'w') as outfile:
        # 写入标题
        outfile.write("sp1\tsp2\tz_value\td_value\n")
        for pair, (z_sum, d_sum, count) in sorted_pairs:
            sp1, sp2 = pair
            # 计算z_value和gamma的平均值
            avg_z = z_sum / count if count > 0 else 0
            avg_d = d_sum / count if count > 0 else 0
            outfile.write(f"{sp1}\t{sp2}\t{avg_z}\t{avg_d}\n")

    print(f"处理完成，结果已保存至 '{output_file}'")

# 从命令行获取输入和输出文件名
if len(sys.argv) != 3:
    print("用法: python species_hybridization.py <输入文件> <输出文件>")
    sys.exit(1)

input_file = sys.argv[1]  # 第一个命令行参数是输入文件
output_file = sys.argv[2]  # 第二个命令行参数是输出文件

# 运行处理函数
process_data(input_file, output_file)
