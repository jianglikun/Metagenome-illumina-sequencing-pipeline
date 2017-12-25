下载地址：
http://soap.genomics.org.cn/soapaligner.html
(安装最新版本2.21，使用2.20的时候报错了)

安装使用：
tar zxvf SOAPaligner.tar.gz
cd SOAPaligner

首先对宿主基因组进行建库：
./soap2.20release/2bwt-builder hg19

首先将paired测序数据合并成为一个：
zcat 170627-1_S6_R* > SBD_input.samples.gz


然后剔除：
../../Screen/soap2.21release/soap -a SBD_input.samples.gz -D ../../Screen/hg19.index -M 4 -l 30 -v 10 -p 8 -o stander_out 2>>log


其中-o得到对文件是比对上宿主的序列
-u得到的序列是未比对上的序列,此数据是我们使用的fasta数据。


因为soap比对后得到的是fasta数据，但是以后的分析如metaphlan2需要fastq数据，此处用python脚本处理比对上的stander_out文件得到，指控后的paired测序数据：

1、python screen_paired_seq.py stander_out 170627-1_S6_R1_001.fastq.trimmed.paired.gz 170627-1_R1.fastq
2、python screen_paired_seq.py stander_out 170627-1_S6_R2_001.fastq.trimmed.paired.gz 170627-1_R2.fastq
