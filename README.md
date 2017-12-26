# 宏基因组分析流程

    宏基因组分析流程来源于Mocat软件，该软件集成了一些二代测序分析软件：如fastqc、metagenemark、cd-hit等等。该流程将肠道菌群的illumina测序数据进行处理：从指控到生成丰度表（物种丰度表、基因丰度表）和基因集。


    Step1：利用fastqc软件查看样本测序结果
    Step2：利用cutaddapter去除序列测序接头
    Step3：利用Solexsaqa去除低质量数据
    Step4：利用soapalign去除宿主reads
    以上为指控过程

    Step5：首先计算组装需要的k-mer
    Step6：利用soapdenove将每个样本的reads进行组装
    Step7：组装得到的scafSeq文件，包含N，有些序列过短，脚本处理得到长度大于500，并且无N的序列
    以上为组装过程

Step8:利用metagenemark将组装结果预测为基因
以上为基因预测

Step9:利用cd-hit生成基因集
以上为生成基因集


Step--run：1、利用指控后的数据与metaphlan2生成物种丰度表
Step--run：2、利用每个样本的基因集注释到功能数据库KEGG、COG、VFDB、ARDB
Step--run：3、利用指控后的数据得到mOTU结果


