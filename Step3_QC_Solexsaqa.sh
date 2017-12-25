下载地址：
http://solexaqa.sourceforge.net/
解压后直接使用

首先剔除低质量reads,质量值大于20
/home/jianglikun/metagenome_pipeline/Linux_x64/SolexaQA++ dynamictrim ../170627-1_S6_R1_001.fastq.gz ../170627-1_S6_R2_001.fastq.gz -h 20 -d ./
然后剔除长度太短的reads，长度大于30
/home/jianglikun/metagenome_pipeline/Linux_x64/SolexaQA++ lengthsort 170627-1_S6_R1_001.fastq.trimmed.gz ./170627-1_S6_R2_001.fastq.trimmed.gz  -l 30 -d ./
