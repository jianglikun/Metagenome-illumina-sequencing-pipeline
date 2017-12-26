虽然我也不清楚为什么要找出motu

网址：http://www.bork.embl.de/software/mOTU/tutorial.motu.standalone.html
可通过网上提供版本下载安装，mOTUs.pl
但是该软件默认进行指控，包括低质量数据等等，使用也很简单，参照官网提供的就可以了。


此处使用同事提供的版本，删除了指控部分：
首先写一个sample文件：
写上样本名就可以了：170306-01

然后，将样本的paired文件连接到指定目录：
if [ "$inpath" != "$readspath" ]; then
   for i in `less $samples`;do mkdir -p $i/reads.processed.solexaqa $i/stats;done &&
   for i in `less $samples`;do ln -s $readspath/$i/reads.processed.solexaqa/$i.pair.1.fq.gz $i/reads.processed.solexaqa; ln -s $readspath/$i/reads.processed.solexaqa/$i.pair.2.fq.gz $i/reads.processed.solexaqa; ln -s $readspath/$i/reads.processed.solexaqa/$i.single.fq.gz $i/reads.processed.solexaqa; ln -s $readspath/$i/stats/$i.readtrimfilter.solexaqa.stats $i/stats;done
fi


--dbdir指定的是motu_data目录：
perl mOTUs.modified.pl --sample-file samples --dbdir /home/jianglikun/metagenome_pipeline/jlk_script/data/test/motus_data/data/ --output-directory ./test --processors 8
