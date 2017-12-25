安装使用详见自己博客：
http://blog.sina.com.cn/s/blog_c4e3e0620102wz87.html

在组装之前要将paired文件cat到一起，并且在config文件中指定。

首先使用soapdenove需要配置config文件：
soap_config_file
配置文件内容：
max_rd_len=150
[LIB]
avg_ins=210
asm_flags=3
rank=1
q1=../screen_170306-01.1.fq.gz
q2=../screen_170306-01.2.fq.gz
[LIB]
asm_flags=1
q=../170306-01.screen.single.fq.gz
此处使用的是soapalign未比对到宿主到fasta序列

/k11e/pvdisk/fastbase/Users/jianglikun/genome/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s ../soap_config_file -k 49 -o soap_runout
