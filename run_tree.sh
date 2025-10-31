# README

## ���
#���������ڴӻ������ƶ������������������ʷ�����̰����������֣�  
#1. **�������ƶ�**������ IQ-TREE ���ɵĻ�������ͨ�� `newick-utils` ���˵� bootstrap ��֧����ʹ�� **ASTRAL** �ƶ���������  
#2. **��״��������**��ʹ�� **PhyloNet** �Բ�ͬ���ݼ����������ƶϣ����ö�� replicate �Ͳ��м�����߽�����Ƚ��ԡ�  

#---

## ��������
#- **IQ-TREE** (������ `.treefile`)  
#- **newick-utils 1.6**  
#- **ASTRAL 5.7.1**  
#- **PhyloNet 3.8.2**  
#- **GNU Parallel 20190622**  
#- **Python 3**��������ű� `generate_input_nex_ml.py`��  
#- **Java 1.8+**  

#---

## ��������

### 1. ����������

# ʹ�� newick-utils ���� IQ-TREE �� bootstrap ֵ
/data/soft/newick-utils/newick-utils-1.6/bin/nw_ed iqtree.trees 'i & b<=33' o > iqtree-BS33.tre

# ʹ�� ASTRAL �ƶ�������
java -jar /data/soft/Astral/ASTRAL-5.7.1/astral.5.7.1.jar -i iqtree-BS33.tre -o astra-BS33.tre

###  2. PhyloNet ��������ƶ�
##(a) Arizonica ���ݼ�

# ��������Ŀ¼
mkdir 00.arizonica
cd 00.arizonica
for i in {1..5}; do mkdir net$i; done

# ���� PhyloNet �����ļ� (5 �� replicate)
for i in {1..5}; do
    python ../generate_input_nex_ml.py arizonica.treefile arizonica_MPL_${i}.nex ${i}
done

# Ϊÿ�� replicate ���� 50 �Σ��� 250 ���ƶ�
for a in {1..5}; do
    for b in {1..50}; do
        echo "java -jar /data/soft/PhyloNet/PhyloNet_3.8.2.jar arizonica_MPL_${a}.nex >net$a/run$b.out" >> run_phylonet.sh
    done
done

# ʹ�� GNU Parallel ����ִ��
/data/soft/parallel/parallel-20190622/src/parallel -j 5 :::: run_phylonet.sh

##(b) Macrogene ���ݼ�
# ��������Ŀ¼
mkdir 01.macrogene
cd 01.macrogene
for i in {1..5}; do mkdir net$i; done

# ���� PhyloNet �����ļ� (5 �� replicate)
for i in {1..5}; do
    python ../generate_input_nex_ml.py macrogene.treefile macrogene_MPL_${i}.nex ${i}
done

# Ϊÿ�� replicate ���� 50 �Σ��� 250 ���ƶ�
for a in {1..5}; do
    for b in {1..50}; do
        echo "java -jar /data/soft/PhyloNet/PhyloNet_3.8.2.jar macrogene_MPL_${a}.nex >net$a/run$b.out" >> run_phylonet.sh
    done
done

# ʹ�� GNU Parallel ����ִ��
/data/soft/parallel/parallel-20190622/src/parallel -j 5 :::: run_phylonet.sh
