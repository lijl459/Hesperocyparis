#!/usr/bin/env python
import sys
from Bio import Phylo
from itertools import combinations

def is_descendant(clade, target):
    """检查 target 是否是 clade 的后代"""
    for descendant in clade.get_terminals():
        if descendant.name == target:
            return True
    return False

def get_species(name):
    """
    提取物种名（去掉 .1, .2 这样的个体编号）
    比如 sp1.1 -> sp1
    """
    return name.split(".")[0]

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <intree> <outname> <outgroup>")
        sys.exit(1)

    intree, outname, outgroup = sys.argv[1], sys.argv[2], sys.argv[3]

    # 读取并解析树文件
    try:
        tree = Phylo.read(intree, "newick")
    except Exception as e:
        print(f"Error reading tree file: {e}")
        sys.exit(1)

    # 获取所有叶子名称
    splist = [term.name for term in tree.get_terminals()]
    lst = []
    if outgroup not in splist:
        sys.exit("{} not in tree".format(outname))

    # 遍历所有不重复的三元组组合
    for sp1, sp2, sp3 in combinations(splist, 3):
        # 判断物种是否不同
        species1, species2, species3 = get_species(sp1), get_species(sp2), get_species(sp3)
        #if len({species1, species2, species3}) < 3:
        #    continue  # 有相同物种（只是不同个体），跳过

        # 三个物种的共同祖先
        anc123 = tree.common_ancestor(sp1, sp2, sp3)

        # 判断外类群是否是 anc123 的后代
        if not is_descendant(anc123, outgroup):
            # 获取两两共同祖先关系，确定排列
            anc12 = tree.common_ancestor(sp1, sp2)
            anc13 = tree.common_ancestor(sp1, sp3)
            anc23 = tree.common_ancestor(sp2, sp3)

            if not is_descendant(anc12, sp3):
                lst.append([sp1, sp2, sp3, outgroup])
            elif not is_descendant(anc13, sp2):
                lst.append([sp1, sp3, sp2, outgroup])
            elif not is_descendant(anc23, sp1):
                lst.append([sp2, sp3, sp1, outgroup])
            # else: 不符合 ((a,b),c) 关系，忽略

    # 将结果写入输出文件
    try:
        with open(outname, "w") as outf:
            for triple in lst:
                outf.write("\t".join(triple) + "\n")
        print(f"Results written to {outname}")
    except Exception as e:
        print(f"Error writing to output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
