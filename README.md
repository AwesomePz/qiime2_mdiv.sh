# qiime2_mdiv
本文档适用于qiime2 amplion基础分析流程：包括降噪、生成特征表、物种注释、树、稀释曲线、a/b多样性

## 参考教程：
- https://amplicon-docs.qiime2.org/en/stable/tutorials/moving-pictures/
- https://amplicon-docs.qiime2.org/en/stable/tutorials/gut-to-soil/

## 使用方式：
1. 上传sample-metadata.csv、treat.csv至$example_dir/00.data（sample-metadata中的样本名必须可以明确指向某一fq）
2. 上传fq.gz至$example_dir/00.data/fq（可分析单双端序列）
3. nohup sh work.sh > qiime2.log 2>&1 &

## 版本：
- qiime2版本：2026.4
- 数据库版本：silva138.2（16S 18S）
