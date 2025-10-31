# README

## 简介
#本流程用于从基因树推断物种树和网络进化历史。流程包含两个部分：  
#1. **物种树推断**：利用 IQ-TREE 生成的基因树，通过 `newick-utils` 过滤低 bootstrap 分支，再使用 **ASTRAL** 推断物种树。  
#2. **网状进化分析**：使用 **PhyloNet** 对不同数据集进行网络推断，采用多次 replicate 和并行计算提高结果的稳健性。  

#---

## 环境需求
#- **IQ-TREE** (已生成 `.treefile`)  
#- **newick-utils 1.6**  
#- **ASTRAL 5.7.1**  
#- **PhyloNet 3.8.2**  
#- **GNU Parallel 20190622**  
#- **Python 3**（需包含脚本 `generate_input_nex_ml.py`）  
#- **Java 1.8+**  

#---

## 分析流程

### 1. 构建物种树

# 使用 newick-utils 过滤 IQ-TREE 的 bootstrap 值
/data/soft/newick-utils/newick-utils-1.6/bin/nw_ed iqtree.trees 'i & b<=33' o > iqtree-BS33.tre

# 使用 ASTRAL 推断物种树
java -jar /data/soft/Astral/ASTRAL-5.7.1/astral.5.7.1.jar -i iqtree-BS33.tre -o astra-BS33.tre

###  2. PhyloNet 网络进化推断
##(a) Arizonica 数据集

# 创建工作目录
mkdir 00.arizonica
cd 00.arizonica
for i in {1..5}; do mkdir net$i; done

# 生成 PhyloNet 输入文件 (5 个 replicate)
for i in {1..5}; do
    python ../generate_input_nex_ml.py arizonica.treefile arizonica_MPL_${i}.nex ${i}
done

# 为每个 replicate 运行 50 次，共 250 次推断
for a in {1..5}; do
    for b in {1..50}; do
        echo "java -jar /data/soft/PhyloNet/PhyloNet_3.8.2.jar arizonica_MPL_${a}.nex >net$a/run$b.out" >> run_phylonet.sh
    done
done

# 使用 GNU Parallel 并行执行
/data/soft/parallel/parallel-20190622/src/parallel -j 5 :::: run_phylonet.sh

##(b) Macrogene 数据集
# 创建工作目录
mkdir 01.macrogene
cd 01.macrogene
for i in {1..5}; do mkdir net$i; done

# 生成 PhyloNet 输入文件 (5 个 replicate)
for i in {1..5}; do
    python ../generate_input_nex_ml.py macrogene.treefile macrogene_MPL_${i}.nex ${i}
done

# 为每个 replicate 运行 50 次，共 250 次推断
for a in {1..5}; do
    for b in {1..50}; do
        echo "java -jar /data/soft/PhyloNet/PhyloNet_3.8.2.jar macrogene_MPL_${a}.nex >net$a/run$b.out" >> run_phylonet.sh
    done
done

# 使用 GNU Parallel 并行执行
/data/soft/parallel/parallel-20190622/src/parallel -j 5 :::: run_phylonet.sh
