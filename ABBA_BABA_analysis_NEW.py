import numpy as np
from Bio import SeqIO
from scipy import stats
import concurrent.futures
import math
import argparse
import csv
import sys
import random

# 碱基编码
BASE_ENCODING = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
MISSING_BASE = 4

def load_fasta(file_path):
    """
    加载FASTA文件，并将序列转换为NumPy数组进行数值编码。
    :param file_path: FASTA文件路径
    :return: 字典，key为物种名称，value为NumPy数组形式的碱基序列
    """
    sequences = {}
    seq_length = None  # 用于检查所有序列长度是否一致
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            # 将碱基转换为大写并编码
            seq = record.seq.upper()
            encoded_seq = np.array([BASE_ENCODING.get(base, MISSING_BASE) for base in seq], dtype=np.int8)
            sequences[record.id] = encoded_seq
            if seq_length is None:
                seq_length = len(encoded_seq)
            elif len(encoded_seq) != seq_length:
                print(f"Error: 序列长度不一致。物种 {record.id} 的序列长度为 {len(encoded_seq)}，预期长度为 {seq_length}.")
                sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    return sequences

def read_species_quadruples(file_path):
    """
    从文件中读取四物种组合，每行包含4个物种名称
    :param file_path: 物种组合文件路径
    :return: 物种四元组列表
    """
    species_quadruples = []
    try:
        with open(file_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                species_quadruple = line.strip().split()
                if len(species_quadruple) == 4:
                    species_quadruples.append(species_quadruple)
                else:
                    print(f"跳过无效的物种组合（行 {line_num}）：{line.strip()}")
    except Exception as e:
        print(f"Error reading species quadruples file: {e}")
        sys.exit(1)
    return species_quadruples

def generate_filtered_sequences(sequences):
    """
    生成一个过滤后的序列，移除包含缺失位点和所有物种碱基相同的位点。
    :param sequences: 字典，key为物种名称，value为NumPy数组形式的碱基序列
    :return: 过滤后的序列字典
    """
    species_list = list(sequences.keys())
    seq_length = len(sequences[species_list[0]])  # 假设所有物种序列长度相同

    # 创建一个二维数组，形状为 (物种数, 序列长度)
    all_seqs = np.vstack([sequences[species] for species in species_list])

    # 检查缺失值
    valid_mask = np.all(all_seqs != MISSING_BASE, axis=0)

    # 检查是否所有物种碱基相同
    same_mask = np.all(all_seqs == all_seqs[0, :], axis=0)

    # 过滤条件：无缺失且不是所有物种碱基相同
    filtered_mask = valid_mask & (~same_mask)

    # 生成新的过滤后的序列
    filtered_sequences = {species: sequences[species][filtered_mask] for species in species_list}

    return filtered_sequences

def calculate_abba_baba_vectorized(sequences, species_quadruple):
    """
    使用向量化方法计算ABBA和BABA计数。
    :param sequences: 字典，key为物种名称，value为NumPy数组形式的碱基序列（已过滤）
    :param species_quadruple: 四物种的组合，例如['species1', 'species2', 'species3', 'species4']
    :return: D值, ABBA计数, BABA计数
    """
    species1, species2, species3, species4 = species_quadruple

    # 获取序列
    try:
        seq1 = sequences[species1]
        seq2 = sequences[species2]
        seq3 = sequences[species3]
        seq4 = sequences[species4]
    except KeyError as e:
        raise ValueError(f"物种 {e} 的序列不存在。")

    # 检查序列长度
    if not (len(seq1) == len(seq2) == len(seq3) == len(seq4)):
        raise ValueError("输入的物种序列长度不一致")

    # ABBA条件：species1 == species4, species2 == species3, species1 != species2
    abba_mask = (seq1 == seq4) & (seq2 == seq3) & (seq1 != seq2)

    # BABA条件：species1 == species3, species2 == species4, species1 != species2
    baba_mask = (seq1 == seq3) & (seq2 == seq4) & (seq1 != seq2)

    # 统计ABBA和BABA
    abba_count = np.sum(abba_mask)
    baba_count = np.sum(baba_mask)

    # 计算D值
    if abba_count + baba_count == 0:
        D = 0.0
    else:
        D = (abba_count - baba_count) / (abba_count + baba_count)

    return D, abba_count, baba_count

def block_jackknife_single(sequences, species_quadruple, block_size, keep_blocks):
    """
    执行单次Block Jackknife重抽样并计算D值。
    :param sequences: 字典，key为物种名称，value为NumPy数组形式的碱基序列
    :param species_quadruple: 四物种组合
    :param block_size: 块大小
    :param keep_blocks: 要保留的块索引列表
    :return: D值
    """
    combined_seq = {}
    for species, seq in sequences.items():
        # 计算每个块的开始和结束位置
        start_positions = [i * block_size for i in keep_blocks]
        end_positions = [start + block_size for start in start_positions]
        # 合并块，确保不超出序列长度
        parts = [seq[start:end] for start, end in zip(start_positions, end_positions) if start < len(seq)]
        if parts:
            combined_seq[species] = np.concatenate(parts)
        else:
            combined_seq[species] = np.array([], dtype=np.int8)

    # 如果任何序列为空，则返回D=0
    if any(len(seq) == 0 for seq in combined_seq.values()):
        return 0.0

    try:
        D, _, _ = calculate_abba_baba_vectorized(combined_seq, species_quadruple)
    except ValueError:
        D = 0.0

    return D

def process_quadruple(sequences, species_quadruple, block_size, num_resamples):
    """
    处理单个物种四元组，计算ABBA、BABA、D值，并进行Block Jackknife重抽样计算Z值和P值。
    :param sequences: 字典，key为物种名称，value为NumPy数组形式的碱基序列
    :param species_quadruple: 四物种组合
    :param block_size: 块大小
    :param num_resamples: 重抽样次数
    :return: 字典包含ABBA计数、BABA计数、D值、Z值、P值
    """
    try:
        # 过滤序列
        seqDict = {species_quadruple[i]: sequences[species_quadruple[i]] for i in range(4)}
        filtered_sequence = generate_filtered_sequences(seqDict)

        # 计算ABBA, BABA, D值
        D_original, ABBA_original, BABA_original = calculate_abba_baba_vectorized(filtered_sequence, species_quadruple)
    except ValueError as e:
        print(f"错误处理物种四元组 {species_quadruple}: {e}")
        return {
            "Species1": species_quadruple[0],
            "Species2": species_quadruple[1],
            "Species3": species_quadruple[2],
            "Species4": species_quadruple[3],
            "ABBA_count": 0,
            "BABA_count": 0,
            "D_value": 0.0,
            "Z_value": 0.0,
            "P_value": 1.0
        }

    if ABBA_original + BABA_original == 0:
        # 如果ABBA和BABA都为0，无法计算D值和后续统计
        return {
            "Species1": species_quadruple[0],
            "Species2": species_quadruple[1],
            "Species3": species_quadruple[2],
            "Species4": species_quadruple[3],
            "ABBA_count": ABBA_original,
            "BABA_count": BABA_original,
            "D_value": 0.0,
            "Z_value": 0.0,
            "P_value": 1.0
        }

    seq_length = len(filtered_sequence[species_quadruple[0]])
    block_count = math.ceil(seq_length / block_size)

    D_resamples = []

    # 进行Block Jackknife重抽样
    for _ in range(num_resamples):
        # 随机选择要保留的块数，这里保留一半到所有块
        keep_block_count = random.randint(math.ceil(int(block_count / 3)), int(block_count * 2 / 3))
        keep_blocks = sorted(random.sample(range(block_count), keep_block_count))
        D = block_jackknife_single(filtered_sequence, species_quadruple, block_size, keep_blocks)
        D_resamples.append(D)

    # 计算标准误差
    if len(D_resamples) > 1:
        SE = np.std(D_resamples, ddof=1) / math.sqrt(len(D_resamples))
    else:
        SE = 0.0

    # 计算Z值
    if SE > 0:
        Z = D_original / SE
    else:
        Z = 0.0

    # 计算P值（双尾检验）
    P = stats.norm.sf(abs(Z)) * 2

    return {
        "Species1": species_quadruple[0],
        "Species2": species_quadruple[1],
        "Species3": species_quadruple[2],
        "Species4": species_quadruple[3],
        "ABBA_count": ABBA_original,
        "BABA_count": BABA_original,
        "D_value": round(D_original, 4),
        "Z_value": round(Z, 4),
        "P_value": round(P, 4)
    }

def block_jackknife_resampling(sequences, species_quadruples, block_size=100, num_resamples=100, num_threads=4):
    """
    执行Block Jackknife重抽样并计算所有D值，使用多进程加速。
    :param sequences: 字典，key为物种名称，value为NumPy数组形式的碱基序列
    :param species_quadruples: 四物种组合的列表
    :param block_size: 块大小
    :param num_resamples: 重抽样次数
    :param num_threads: 进程数量
    :return: 列表，每个元素是一个字典包含一个物种四元组的结果
    """
    results = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        # 提交所有任务
        futures = [
            executor.submit(
                process_quadruple,
                sequences,
                quad,
                block_size,
                num_resamples, 
            )
            for quad in species_quadruples
        ]

        # 收集结果
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.append(result)

    return results

def save_results(results, output_file):
    """
    将计算结果保存到CSV文件中
    :param results: 结果列表
    :param output_file: 输出文件路径
    """
    try:
        with open(output_file, mode="w", newline="") as file:
            fieldnames = ["Species1", "Species2", "Species3", "Species4",
                          "ABBA_count", "BABA_count", "D_value",
                          "Z_value", "P_value"]
            writer = csv.DictWriter(file, fieldnames=fieldnames)

            writer.writeheader()
            for row in results:
                writer.writerow(row)
    except Exception as e:
        print(f"Error saving results to {output_file}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="计算ABBA, BABA, D值及Block Jackknife分析")
    parser.add_argument("-f", "--fasta", required=True, help="输入FASTA文件路径")
    parser.add_argument("-s", "--species_quadruples", required=True, help="物种四元组文件路径")
    parser.add_argument("-b", "--block_size", type=int, default=100, help="块大小，默认值为100")
    parser.add_argument("-r", "--resamples", type=int, default=100, help="重抽样次数，默认值为100")
    parser.add_argument("-t", "--threads", type=int, default=4, help="进程数量，默认值为4")
    parser.add_argument("-o", "--output", required=True, help="输出结果文件路径（CSV格式）")
    args = parser.parse_args()

    fasta_file = args.fasta
    species_quadruples_file = args.species_quadruples
    block_size = args.block_size
    num_resamples = args.resamples
    num_threads = args.threads
    output_file = args.output

    print("加载FASTA序列...")
    sequences = load_fasta(fasta_file)
    print(f"加载完成，共加载到 {len(sequences)} 个物种序列。")

    print("读取物种四元组...")
    species_quadruples = read_species_quadruples(species_quadruples_file)
    print(f"读取完成，共有 {len(species_quadruples)} 个物种四元组。")

    print(f"开始进行Block Jackknife重抽样（块大小={block_size}, 重抽样次数={num_resamples}, 进程数={num_threads})...")
    results = block_jackknife_resampling(
        sequences,
        species_quadruples,
        block_size=block_size,
        num_resamples=num_resamples,  # 修正后的关键字参数
        num_threads=num_threads
    )
    print("Block Jackknife重抽样完成。")

    print("准备结果并写入文件...")
    # 写入CSV文件
    if results:
        try:
            save_results(results, output_file)
            print(f"结果文件写入完成：{output_file}")
        except Exception as e:
            print(f"Error writing output file: {e}")
            sys.exit(1)
    else:
        print("没有结果需要写入。")

    print("所有操作完成。")

if __name__ == "__main__":
    # 设置随机种子以确保可重复性
    random.seed(42)
    np.random.seed(42)
    main()
