下载安装：
git clone https://github.com/marcelm/cutadapt/
          cd cutadapt
          python setup.py build_ext -i
          cd bin
          export PATH=$PATH:$PWD

使用：
Adapter=`grep -v ">" ../adapter.fa|head -1`
AdapterRead2=`grep -v ">" ../adapter.fa|tail -1`
注意检查adapter变量是否成功建立。

cutadapt -a $Adapter -A $AdapterRead2 -m 30 -o ./170627-1.trimadapter.1.fq.gz -p ./170627-1.trimadapter.2.fq.gz ./170627-1.1.fq.gz ./170627-1.2.fq.gz > 170627-1.cutadapt.stats
