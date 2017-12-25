下载
    $ wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip
解压
    $ unzip fastqc_v0.11.3.zip
设置权限
    $ cd FastQC/
    $ chmod 755 fastqc
加入到 PATH
    $ export PATH=/home/user/FastQC/:$PATH
测试
    $ fastqc --help

使用:
1、首先将两个paired文件cat到一起
2、fastqc -o $sample --noextract -t 8 $sample.raw.fq >>log
