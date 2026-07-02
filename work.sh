# 微生物多样性流程分析（测试）

# 准备文件
# 00.data/fq：存放所有fq.gz
# 00.data/sample-metadata.csv：样品编号	分组
# 00.data/treat.csv：要比对的分组

set -e
source /home/zpy/software/miniconda3/etc/profile.d/conda.sh
conda activate qiime2

# 00 输入文件
WORK_DIR=/home/zpy/Project/mdiv/test_MbPL202511786
cd $WORK_DIR
sed 's/\r//g; s/,/\t/g' "$WORK_DIR/00.data/sample-metadata.csv" > "$WORK_DIR/00.data/sample-metadata.tsv"
sed -i '1i sample-id\tGroup\n#q2:types\tcategorical' "$WORK_DIR/00.data/sample-metadata.tsv"
sed -i -e '$a\' "$WORK_DIR/00.data/sample-metadata.tsv"
sed 's/\r//g; s/,/\t/g' "$WORK_DIR/00.data/treat.csv" > "$WORK_DIR/00.data/treat.tsv"

# fastp
FQ_DIR="$WORK_DIR/00.data/fq"
CLEAN_DIR="$WORK_DIR/00.data/fq_clean"
mkdir -p "$CLEAN_DIR"
for r1 in "$FQ_DIR"/*_R1.fq.gz; do
    [ -f "$r1" ] || continue
    sample=$(basename "$r1" | sed 's/_R1\.fq\.gz//')
    r2="$FQ_DIR/${sample}_R2.fq.gz"
    if [ -f "$r2" ]; then
        fastp -i "$r1" -I "$r2" -o "$CLEAN_DIR/${sample}_R1.clean.fq.gz"  -O "$CLEAN_DIR/${sample}_R2.clean.fq.gz" -q 20 -u 30 -l 50 --thread 4
    else
        fastp -i "$r1" -o "$CLEAN_DIR/${sample}_R1.clean.fq.gz" -q 20 -u 30 -l 50 --thread 4
    fi
done
mv fastp.* $CLEAN_DIR

# 生成sample-metadata-treat*.tsv
META=$WORK_DIR/00.data/sample-metadata.tsv
TREAT=$WORK_DIR/00.data/treat.tsv
t=1
while IFS=$'\t' read -a cols; do
    [ -z "${cols[0]}" ] && continue
    f="$WORK_DIR/00.data/sample-metadata-treat${t}.tsv"
    echo -e "sample-id\tGroup\n#q2:types\tcategorical" > "$f"
    while IFS=$'\t' read sid g; do
        for c in "${cols[@]}"; do
            [ "$g" == "$c" ] && echo -e "$sid\t$g" && break
        done
    done < "$META" >> "$f"
    ((t++))
done < "$TREAT"

# 生成fq-manifest.tsv
MANIFEST="$WORK_DIR/00.data/fq-manifest.tsv"
CLEAN_DIR="$WORK_DIR/00.data/fq_clean"
mapfile -t meta_ids < <(awk -F'\t' 'NR>2 {print $1}' "$META")
has_r2=$(ls "$CLEAN_DIR"/*_R2.clean.fq.gz 2>/dev/null | head -1)
if [ -n "$has_r2" ]; then
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"
else
    echo -e "sample-id\tabsolute-filepath" > "$MANIFEST"
fi
for sid in "${meta_ids[@]}"; do
    found=0
    for r1_file in "$CLEAN_DIR"/*_R1.clean.fq.gz; do
        [ -f "$r1_file" ] || continue
        basename=$(basename "$r1_file" _R1.clean.fq.gz)
        clean_sid="${basename##*-}"
        if [ "$clean_sid" == "$sid" ]; then
            r1=$(realpath "$r1_file")    
            if [ -n "$has_r2" ]; then
                r2_file="${r1_file/_R1.clean.fq.gz/_R2.clean.fq.gz}"
                if [ -f "$r2_file" ]; then
                    r2=$(realpath "$r2_file")
                    echo -e "$sid\t$r1\t$r2" >> "$MANIFEST"
                    found=1
                    break
                fi
            else
                echo -e "$sid\t$r1" >> "$MANIFEST"
                found=1
                break
            fi
        fi
    done
    if [ "$found" -eq 0 ]; then
        echo "警告: 未找到样本 $sid 的 clean 文件" >&2
    fi
done

cols=$(head -1 "$MANIFEST" | awk -F'\t' '{print NF}')
if [ "$cols" -eq 3 ]; then
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST" \
        --output-path 00.data/demux.qza \
        --input-format PairedEndFastqManifestPhred33V2
else
    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path "$MANIFEST" \
        --output-path 00.data/demux.qza \
        --input-format SingleEndFastqManifestPhred33V2
fi

qiime demux summarize \
  --i-data 00.data/demux.qza \
  --o-visualization 00.data/demux.qzv

# 01 质控+降噪+拼接
mkdir 01.dada
if [ "$cols" -eq 3 ]; then
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs 00.data/demux.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads 0 \
        --o-representative-sequences 01.dada/asv-seqs.qza \
        --o-table 01.dada/asv-table.qza \
        --o-denoising-stats 01.dada/denoising-stats.qza \
        --o-base-transition-stats 01.dada/base-transition-stats.qza
else
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs 00.data/demux.qza \
        --p-trunc-len 0 \
        --p-n-threads 0 \
        --o-representative-sequences 01.dada/asv-seqs.qza \
        --o-table 01.dada/asv-table.qza \
        --o-denoising-stats 01.dada/denoising-stats.qza \
        --o-base-transition-stats 01.dada/base-transition-stats.qza
fi #耗时

# 重命名
mkdir -p 01.dada/renamed
qiime tools export \
  --input-path 01.dada/asv-seqs.qza \
  --output-path 01.dada/renamed
awk '/^>/{printf ">ASV_%d\n", ++n; next} {print}' 01.dada/renamed/dna-sequences.fasta > 01.dada/renamed/asv_renamed.fasta
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 01.dada/renamed/asv_renamed.fasta \
  --output-path 01.dada/asv-seqs.qza
qiime tools export \
  --input-path 01.dada/asv-table.qza \
  --output-path 01.dada/renamed
biom convert \
  -i 01.dada/renamed/feature-table.biom \
  -o 01.dada/renamed/feature-table.tsv \
  --to-tsv
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 && $0 ~ /^#/ && $0 !~ /^#OTU ID/ {print; next} /^#OTU ID/ {gsub(/^#OTU ID/,"ASV_ID",$1); print $0; next} {old=$1; $1="ASV_" ++n; map[old]=$1; print $0} END {print "ASV_ID\tsequence" > "01.dada/asv-id-map.tsv"; for(old in map) print map[old] "\t" old > "01.dada/asv-id-map.tsv"}' 01.dada/renamed/feature-table.tsv > 01.dada/renamed/feature-table-renamed.tsv
biom convert \
  -i 01.dada/renamed/feature-table-renamed.tsv \
  -o 01.dada/renamed/feature-table-renamed.biom \
  --to-hdf5
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path 01.dada/renamed/feature-table-renamed.biom \
  --output-path 01.dada/asv-table.qza
rm -rf 01.dada/renamed

# 02 生成特征表和代表序列
mkdir 02.feature
qiime feature-table filter-features \
  --i-table 01.dada/asv-table.qza \
  --p-min-frequency 2 \
  --o-filtered-table 02.feature/asv-table-mf2.qza
qiime feature-table filter-seqs \
  --i-data 01.dada/asv-seqs.qza \
  --i-table 02.feature/asv-table-mf2.qza \
  --o-filtered-data 02.feature/asv-seqs-mf2.qza
qiime feature-table summarize \
  --i-table 02.feature/asv-table-mf2.qza \
  --m-metadata-file 00.data/sample-metadata.tsv \
  --o-summary 02.feature/asv-table-mf2.qzv \
  --o-sample-frequencies 02.feature/sample-frequencies-mf2.qza \
  --o-feature-frequencies 02.feature/asv-frequencies-mf2.qza

for meta_file in 00.data/sample-metadata-treat*.tsv; do
    treat=$(basename "$meta_file" | sed 's/sample-metadata-treat//; s/\.tsv$//')
    mkdir -p "02.feature/treat${treat}/"
    qiime feature-table filter-samples \
        --i-table 01.dada/asv-table.qza \
        --m-metadata-file "$meta_file" \
        --o-filtered-table "02.feature/treat${treat}/asv-table-treat${treat}.qza"
    qiime feature-table filter-features \
        --i-table "02.feature/treat${treat}/asv-table-treat${treat}.qza" \
        --p-min-frequency 2 \
        --o-filtered-table "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qza"
    qiime feature-table filter-seqs \
        --i-data 01.dada/asv-seqs.qza \
        --i-table "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qza" \
        --o-filtered-data "02.feature/treat${treat}/asv-seqs-mf2-treat${treat}.qza"
    qiime feature-table summarize \
        --i-table "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qza" \
        --m-metadata-file "$meta_file" \
        --o-summary "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qzv" \
        --o-sample-frequencies "02.feature/treat${treat}/sample-frequencies-mf2-treat${treat}.qza" \
        --o-feature-frequencies "02.feature/treat${treat}/asv-frequencies-mf2-treat${treat}.qza"
done

# 03 物种分类
mkdir 03.taxonomy
CLASSIFIER="/home/zpy/software/qiime2/silva138.2"
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads 02.feature/asv-seqs-mf2.qza \
  --o-classification 03.taxonomy/taxonomy.qza #耗时
qiime feature-table tabulate-seqs \
  --i-data 02.feature/asv-seqs-mf2.qza \
  --i-taxonomy 03.taxonomy/taxonomy.qza \
  --m-metadata-file 02.feature/asv-frequencies-mf2.qza \
  --o-visualization 03.taxonomy/asv-seqs-mf2.qzv
qiime taxa barplot \
  --i-table 02.feature/asv-table-mf2.qza \
  --i-taxonomy 03.taxonomy/taxonomy.qza \
  --m-metadata-file 00.data/sample-metadata.tsv \
  --o-visualization 03.taxonomy/taxa-bar-plots.qzv

# ASV丰度表
qiime tools export \
  --input-path 02.feature/asv-table-mf2.qza \
  --output-path /tmp/asv-table-mf2
biom convert \
  -i /tmp/asv-table-mf2/feature-table.biom \
  -o /tmp/asv-table-mf2/feature-table.tsv \
  --to-tsv
qiime tools export \
  --input-path 03.taxonomy/taxonomy.qza \
  --output-path /tmp/tax
awk -F '\t' 'BEGIN{OFS="\t"} NR==FNR && FNR>1 {tax[$1]=$2; next} FNR==1 {next} FNR==2 {gsub(/^#OTU ID/,"ASV_ID",$1); print $0,"taxonomy"; next}{print $0,tax[$1]}' /tmp/tax/taxonomy.tsv /tmp/asv-table-mf2/feature-table.tsv > 02.feature/feature-tax.tsv
rm -rf /tmp/tax

# 04 树
mkdir 04.tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 02.feature/asv-seqs-mf2.qza \
  --p-n-threads 0 \
  --o-alignment 04.tree/aligned-rep-seqs.qza \
  --o-masked-alignment 04.tree/masked-aligned-rep-seqs.qza \
  --o-tree 04.tree/unrooted-tree.qza \
  --o-rooted-tree 04.tree/rooted-tree.qza #耗时

# 05 稀释曲线
mkdir 05.alpha-rarefaction
MEDIAN_DEPTH=$(awk -F'\t' 'NR>2 {for(i=2; i<=NF; i++) sum[i]+=$i} END {n=0; for(i=2; i<=NF; i++) arr[++n]=sum[i]; asort(arr); if(n%2) print arr[int(n/2)+1]; else print int((arr[n/2]+arr[n/2+1])/2)}' /tmp/asv-table-mf2/feature-table.tsv)
qiime diversity alpha-rarefaction \
  --i-table 02.feature/asv-table-mf2.qza \
  --i-phylogeny 04.tree/rooted-tree.qza \
  --p-min-depth 10 \
  --p-max-depth $MEDIAN_DEPTH \
  --m-metadata-file 00.data/sample-metadata.tsv \
  --o-visualization 05.alpha-rarefaction/alpha-rarefaction.qzv #--p-max-depth依据02.feature/asv-table-ms2.qza中50%中位数确定
rm -rf /tmp/asv-table-mf2/

for meta_file in 00.data/sample-metadata-treat*.tsv; do
    treat=$(basename "$meta_file" | sed 's/sample-metadata-treat//; s/\.tsv$//')
    qiime tools export \
      --input-path "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qza" \
      --output-path "/tmp/asv-table-mf2-treat${treat}"
    biom convert \
      -i "/tmp/asv-table-mf2-treat${treat}/feature-table.biom" \
      -o "/tmp/asv-table-mf2-treat${treat}/feature-table.tsv" \
      --to-tsv
    MEDIAN_DEPTH=$(awk -F'\t' 'NR>2 {for(i=2; i<=NF; i++) sum[i]+=$i} END {n=0; for(i=2; i<=NF; i++) arr[++n]=sum[i]; asort(arr); if(n%2) print arr[int(n/2)+1]; else print int((arr[n/2]+arr[n/2+1])/2)}' "/tmp/asv-table-mf2-treat${treat}/feature-table.tsv")
    echo "treat${treat} median depth: $MEDIAN_DEPTH"
    qiime diversity alpha-rarefaction \
      --i-table "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qza" \
      --i-phylogeny 04.tree/rooted-tree.qza \
      --p-min-depth 10 \
      --p-max-depth "$MEDIAN_DEPTH" \
      --m-metadata-file "$meta_file" \
      --o-visualization "05.alpha-rarefaction/alpha-rarefaction-treat${treat}.qzv"
done

# 06 多样性
mkdir -p 06.diversity
for meta_file in 00.data/sample-metadata-treat*.tsv; do
    treat=$(basename "$meta_file" | sed 's/sample-metadata-treat//; s/\.tsv$//')
    RARE_DEPTH=$(awk -F'\t' 'NR>2 {for(i=2;i<=NF;i++) sum[i]+=$i} END {min=sum[2]; for(i=2;i<=NF;i++) if(sum[i]<min) min=sum[i]; print int(min*0.95)+1}' "/tmp/asv-table-mf2-treat${treat}/feature-table.tsv") # 计算rarefaction深度（最低样本深度的95%，向上取整）
    echo "treat${treat} rare depth: $RARE_DEPTH"

    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny 04.tree/rooted-tree.qza \
      --i-table "02.feature/treat${treat}/asv-table-mf2-treat${treat}.qza" \
      --p-sampling-depth "$RARE_DEPTH" \
      --m-metadata-file "$meta_file" \
      --output-dir "06.diversity/treat${treat}/1.divdata"
    rm -rf /tmp/asv-table-mf2-treat${treat}
done

for meta_file in 00.data/sample-metadata-treat*.tsv; do
    treat=$(basename "$meta_file" | sed 's/sample-metadata-treat//; s/\.tsv$//')
    mkdir -p "06.diversity/treat${treat}/2.a-div"
    qiime diversity alpha-group-significance \
      --i-alpha-diversity "06.diversity/treat${treat}/1.divdata/faith_pd_vector.qza" \
      --m-metadata-file "$meta_file" \
      --o-visualization "06.diversity/treat${treat}/2.a-div/faith-pd-group-significance.qzv"
    qiime diversity alpha-group-significance \
      --i-alpha-diversity "06.diversity/treat${treat}/1.divdata/evenness_vector.qza" \
      --m-metadata-file "$meta_file" \
      --o-visualization "06.diversity/treat${treat}/2.a-div/evenness-group-significance.qzv"
    qiime diversity alpha-group-significance \
      --i-alpha-diversity "06.diversity/treat${treat}/1.divdata/shannon_vector.qza" \
      --m-metadata-file "$meta_file" \
      --o-visualization "06.diversity/treat${treat}/2.a-div/shannon-group-significance.qzv"    
    qiime diversity alpha-group-significance \
      --i-alpha-diversity "06.diversity/treat${treat}/1.divdata/observed_features_vector.qza" \
      --m-metadata-file "$meta_file" \
      --o-visualization "06.diversity/treat${treat}/2.a-div/observed-features-group-significance.qzv"

    mkdir -p "06.diversity/treat${treat}/3.b-div"
    qiime diversity beta-group-significance \
      --i-distance-matrix "06.diversity/treat${treat}/1.divdata/unweighted_unifrac_distance_matrix.qza" \
      --m-metadata-file "$meta_file" \
      --m-metadata-column Group \
      --p-pairwise \
      --o-visualization "06.diversity/treat${treat}/3.b-div/unweighted-unifrac-group-significance.qzv" #--m-metadata-column为分组信息，属于00.data/sample-metadata.tsv某列；--o-visualization的文件名记得修改为对应分组
    qiime diversity beta-group-significance \
      --i-distance-matrix "06.diversity/treat${treat}/1.divdata/weighted_unifrac_distance_matrix.qza" \
      --m-metadata-file "$meta_file" \
      --m-metadata-column Group \
      --p-pairwise \
      --o-visualization "06.diversity/treat${treat}/3.b-div/weighted-unifrac-group-significance.qzv"
    qiime diversity beta-group-significance \
      --i-distance-matrix "06.diversity/treat${treat}/1.divdata/bray_curtis_distance_matrix.qza" \
      --m-metadata-file "$meta_file" \
      --m-metadata-column Group \
      --p-pairwise \
      --o-visualization "06.diversity/treat${treat}/3.b-div/bray-curtis-group-significance.qzv"
    qiime diversity beta-group-significance \
      --i-distance-matrix "06.diversity/treat${treat}/1.divdata/jaccard_distance_matrix.qza" \
      --m-metadata-file "$meta_file" \
      --m-metadata-column Group \
      --p-pairwise \
      --o-visualization "06.diversity/treat${treat}/3.b-div/jaccard-group-significance.qzv"
    rm -rf "/tmp/asv-table-mf2-treat${treat}"
done

conda deactivate
