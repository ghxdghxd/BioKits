# install

```
pip install qiime
print_qiime_config.py -t
wget <https://github.com/qiime/qiime-deploy/archive/master.zip>
unzip master.zip
mv qiime-deploy-master qiime-deploy
wget <https://github.com/qiime/qiime-deploy-conf/archive/master.zip>
unzip master.zip
mv qiime-deploy-conf-master qiime-deploy-conf
cd qiime-depoly
python qiime-depoly.py $PATH/qiime_software -f qiime-deploy-
conf/qiime-1.9.1/qiime.conf \--force-remove-failed-dirs
```
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
source activate qiime1 激活qiime1环境
print_qiime_config.py -t
source deactivate 退出qiime1环境
conda remove --name qiime1 --all 删除qiime1环境

# fastq-join

  File "/share/apps/anaconda/bin/print_qiime_config.py", line 40, in
&lt;module&gt;

    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

ImportError: No module named externals

    conda install -c https://conda.anaconda.org/biocore scikit-bio

## error

Additional information: on this system no suitable UTF-8
locales were discovered.  This most likely requires resolving
by reconfiguring the locale system.

解决 export LC_ALL=en_US.utf8