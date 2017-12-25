将指控之后的测序数据进行组装之前需要先计算Kmer
需要输入文件为压缩之后的。

利用现有的perl脚本:详见文章后部。

perl kmer.count.pl ./data/170627-1_S6_R1_001.fastq.trimmed.paired.gz ./data/170627-1_S6_R2_001.fastq.trimmed.paired.gz

直接输出结果：49
