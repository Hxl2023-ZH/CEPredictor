#!/usr/bin/env bash

set -euo pipefail

############################################
# 默认参数
############################################
PROMOTER_LEN=300

############################################
# 帮助信息
############################################
usage() {
    echo "Usage: $0 -g <gff_dir> -f <fna_dir> -r <rbox_table> [-o <out_dir>] [-l <promoter_length>]"
    echo ""
    echo "Required:"
    echo "  -g   GFF文件目录"
    echo "  -f   基因组fna目录"
    echo "  -r   rBOX注释表（包含基因ID）"
    echo ""
    echo "Optional:"
    echo "  -o   输出目录 (default: ./rbox_promoter_result)"
    echo "  -l   promoter长度 (default: 300)"
    exit 1
}

############################################
# 参数解析
GFF_DIR=""
FNA_DIR=""
RBOX_TABLE=""
OUT_DIR="./rbox_promoter_result"

while getopts "g:f:r:o:l:h" opt; do
    case $opt in
        g) GFF_DIR=$OPTARG ;;
        f) FNA_DIR=$OPTARG ;;
        r) RBOX_TABLE=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        l) PROMOTER_LEN=$OPTARG ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$GFF_DIR" || -z "$FNA_DIR" || -z "$RBOX_TABLE" ]]; then
    usage
fi

############################################
# 目录结构
mkdir -p ${OUT_DIR}/{bed,promoter_fasta,background,rbox_promoter,tmp}

############################################
# Step 0: 检查依赖
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found"; exit 1; }
command -v seqkit >/dev/null 2>&1 || { echo "seqkit not found"; exit 1; }
command -v fasta-get-markov >/dev/null 2>&1 || { echo "MEME suite not found"; exit 1; }

############################################
# Step 1: 处理每个GFF
for gff in ${GFF_DIR}/*.gff; do
    [ -e "$gff" ] || continue

    base=$(basename "$gff")
    prefix=${base%.gff}

    echo "===== Processing $prefix ====="

    genome=${FNA_DIR}/${prefix}.fna

    if [[ ! -f "$genome" ]]; then
        echo "[WARNING] Genome not found: $genome"
        continue
    fi

    # Step 1.1: 建立fai索引
    if [[ ! -f "${genome}.fai" ]]; then
        echo "[INFO] Indexing genome..."
        samtools faidx "$genome"
    fi

    # Step 1.2: GFF → BED（CDS）
    bed=${OUT_DIR}/bed/${prefix}.bed

    awk '$3=="CDS"' "$gff" | \
    awk 'BEGIN{OFS="\t"} {
        split($9,a,";");
        id=a[1];
        sub("ID=","",id);
        print $1,$4-1,$5,id,0,$7
    }' > "$bed"

    # Step 1.3: 提取promoter
    promoter_bed=${OUT_DIR}/bed/${prefix}_promoter.bed

    bedtools flank \
        -i "$bed" \
        -g "${genome}.fai" \
        -l ${PROMOTER_LEN} -r 0 \
        > "$promoter_bed"

    # Step 1.4: 提取序列
    promoter_fa=${OUT_DIR}/promoter_fasta/${prefix}_promoter.fasta

    bedtools getfasta \
        -fi "$genome" \
        -bed "$promoter_bed" \
        -s -name \
        > "$promoter_fa"

    # Step 1.5: 背景模型
    fasta-get-markov \
        -dna "$promoter_fa" \
        ${OUT_DIR}/background/${prefix}_background.txt

    # Step 1.6: 提取rBOX对应promoter

    # 从rBOX表中筛选当前MAG的基因ID
    rbox_ids=${OUT_DIR}/tmp/${prefix}_rbox_id.txt

    awk -v p="$prefix" '$0 ~ p {print $1}' "$RBOX_TABLE" > "$rbox_ids"

    if [[ ! -s "$rbox_ids" ]]; then
        echo "[WARNING] No rBOX genes found in $prefix"
        continue
    fi

    # 提取对应序列
    seqkit grep \
        -n -r -f "$rbox_ids" \
        "$promoter_fa" \
        > ${OUT_DIR}/rbox_promoter/${prefix}_rbox.fasta

    echo "[DONE] $prefix"

done

echo "=== All finished ==="
