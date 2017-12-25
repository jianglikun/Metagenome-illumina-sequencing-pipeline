下载安装：
下载地址http://exon.gatech.edu/Genemark/license_download.cgi

除下载软件之外还有一个gm_key_64

1、直接解压MetaGeneMark_linux_64.tar.gz

2、解压gm_key_64，并且  cp gm_key ~/.gm_key

之后就可以直接使用：
~/metagenome_pipeline/jlk_script/metagenemark/MetaGeneMark_linux_64/mgm/gmhmmp -a -d -f L -m ~/metagenome_pipeline/jlk_script/metagenemark/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod -o ./metagenemark.list S030.scaftig -A metagenemark_pro -D metagenemark_nul

蛋白质序列：metagenemark_pro
核算序列：metagenemark_nul
