#! /usr/bin/bash

################# Information ###############
# Name : Qiaozhenghao
# Data : 2024-12-23
# Email: qzhrookie@163.com
#############################################

getopt_cmd=$(getopt -o m:g:2:1:o:h --long help -n $(basename $0) -- "$@")
[ $? -ne 0 ] && exit 1
eval set -- "$getopt_cmd"
# 解析选项
while [ -n "$1" ]
do
    case "$1" in
        -1)
            fq1="$2"
            shift ;;
        -2)
            fq2="$2"
            shift ;;
        -g)
            species="$2"
            shift ;;
        -o)
            out_dir="$2"
            shift ;;
        -m)
            method="$2"
            shift ;;
		-h|--help)
            echo -e "$help_str"
            exit ;;
        --) shift
            break ;;
         *) echo "$1 is not an option"
            exit 1 ;;
    esac
    shift
done

if [ "$species" == mouse ]
then
    genome="/home/qzh/ref/mm10/hisat2_index/ensemble/genome"
    rRNA_index="/home/qzh/ref/mm10/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/mm10/fasta/mm10_ensemble.gtf"
    ID_trans_file="/home/qzh/ref/mm10/Ensemble_id_trans.txt"
    TEgtf="/home/qzh/ref/mm10/TEtranscript_ref/mm10_rmsk_TE.gtf"
elif [ "$species" == pig ]
then
    genome="/home/qzh/ref/ss11/hisat2_index/genome/genome"
    rRNA_index="/home/qzh/ref/ss11/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/ss11/fasta/ss11_genome.gtf"
    ID_trans_file="/home/qzh/ref/ss11/Ensemble_id_trans.txt"
elif [ "$species" == cow ]
then
	genome="/home/qzh/ref/ss11/hisat2_index/genome/genome"
	rRNA_index="/home/qzh/ref/ss11/hisat2_index/rRNA/rRNA"
	gtf_file="/home/qzh/ref/ss11/fasta/ss11_genome.gtf"
	ID_trans_file="/home/qzh/ref/ss11/Ensemble_id_trans.txt"
elif [ "$species" == rat ]
then
    genome="/home/qzh/ref/rn7/hisat2_index/genome/genome"
    rRNA_index="/home/qzh/ref/rn7/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/rn7/genome_annotation_data/rn7_ensemble.gtf"
    ID_trans_file="/home/qzh/ref/rn7/genome_annotation_data/Ensemble_id_trans.txt"
elif [ "$species" == macaca ]
then
    genome="/home/qzh/ref/Macaca/hisat2_index/genome/genome"
    rRNA_index="/home/qzh/ref/Macaca/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/Macaca/fasta/Macaca_ensemble.gtf"
    ID_trans_file="/home/qzh/ref/Macaca/annotation/Ensemble_id_trans.txt"
elif [ "$species" == rabbit ]
then
    genome="/home/qzh/ref/OryCun2/hisat2_index/genome/genome"
    rRNA_index="/home/qzh/ref/OryCun2/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/OryCun2/fasta/OryCun2_genome_ensemble.gtf"
    ID_trans_file="/home/qzh/ref/OryCun2/fasta/Ensemble_id_trans.txt"
elif [ "$species" == human ]
then
    genome="/home/qzh/ref/hg38/hisat2_index/genome/genome"
    rRNA_index="/home/qzh/ref/hg38/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/hg38/fasta/hg38_genome_ensemble.gtf"
    ID_trans_file="/home/qzh/ref/hg38/genome_annotation_data/Ensemble_id_trans.txt"
elif [ "$species" == hamster ]
then
    genome="/home/qzh/ref/MesAur1.0/hisat2_index/genome/genome"
    rRNA_index="/home/qzh/ref/MesAur1.0/hisat2_index/rRNA/rRNA"
    gtf_file="/home/qzh/ref/MesAur1.0/fasta/MesAur_ensemble_genome.gtf"
    ID_trans_file="/home/qzh/ref/MesAur1.0/genome_annotation_data/Ensemble_id_trans.txt"
else
    exit 1;
fi

function Smart3_clean() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

    local fq1_clean="${out_dir}/${fname}_1_val_1.fq.gz"
    local fq2_clean="${out_dir}/${fname}_2_val_2.fq.gz"
    local C_clean_log="${out_dir}/${fname}_C.log"

    local fq1_A_clean="${out_dir}/${fname}_A_1.fq.gz"
    local fq2_A_clean="${out_dir}/${fname}_A_2.fq.gz"
    local A_clean_log="${out_dir}/${fname}_A.log"
    [[ -f ${fq1_clean} ]] || \
        cutadapt -j 20 -n 2 -m 30 -e 0.2 \
            -g AGATGTGTATAAGAGACAG \
            -G AGATGTGTATAAGAGACAG \
            -o ${fq1_A_clean}  \
            -p ${fq2_A_clean} \
            ${fq1} ${fq2} \
            1> ${A_clean_log}

    local fq1_B_clean="${out_dir}/${fname}_B_1.fq.gz"
    local fq2_B_clean="${out_dir}/${fname}_B_2.fq.gz"
    local B_clean_log="${out_dir}/${fname}_B.log"        
    [[ -f ${fq1_clean} ]] || \
        cutadapt -j 20 -n 2 -m 30 -e 0.2 \
            -a CTGTCTCTTATACACATCT \
            -A CTGTCTCTTATACACATCT \
            -o ${fq1_B_clean}  \
            -p ${fq2_B_clean} \
            ${fq1_A_clean} ${fq2_A_clean} \
            1> ${B_clean_log}

    [[ -f ${fq1_clean} ]] || \
        cutadapt -j 20 -n 2 -m 30 -O 5 -e 0.2 \
            -b AAGCAGTGGTATCAACGCAGAGTAC \
            -B AAGCAGTGGTATCAACGCAGAGTAC \
            -a A{50} -a C{50} -a G{50} -a T{50} \
            -A A{50} -A C{50} -A G{50} -A T{50} \
            -o ${fq1_clean} \
            -p ${fq2_clean} \
            ${fq1_B_clean} ${fq2_B_clean} \
            > ${C_clean_log}

    wait
    # Remove
    [[ -f ${fq1_B_clean} ]] && rm ${fq1_B_clean}
    [[ -f ${fq2_B_clean} ]] && rm ${fq2_B_clean}
    [[ -f ${fq1_A_clean} ]] && rm ${fq1_A_clean}
    [[ -f ${fq2_A_clean} ]] && rm ${fq2_A_clean}    

}
function Smart2_clean() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

    local fq1_clean2="${out_dir}/${fname}_1_val_1.fq.gz"
    [[ -f ${fq1_clean2} ]] || \
        trim_galore -q 20 \
            --length 30 -j 8 \
			-a A{15} -a T{15} -a G{15} -a C{15} \
			-a2 A{15} -a2 T{15} -a2 G{15} -a2 C{15}  \
            --paired --clip_R1 10 --clip_R2 10 \
            --three_prime_clip_R1 3 --three_prime_clip_R2 3 \
            ${fq1} ${fq2} \
            -o ${out_dir}
}
function align_genome() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1_*})

    mkdir -p ${out_dir}/rRNA
    mkdir -p ${out_dir}/genome
    local bam_q10="${out_dir}/genome/${fname}_genome_Q10.bam"

    # rRNA
    local rRNA_sam="${out_dir}/rRNA/${fname}_rRNA.sam"
    local rRNA_unmap="${out_dir}/rRNA/${fname}_rRNA_unmap"
    local rRNA_log="${out_dir}/rRNA/${fname}_rRNA.log"
    [[ -f ${bam_q10} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${rRNA_index} \
            --un-conc-gz ${rRNA_unmap} \
            -1 ${fq1} \
            -2 ${fq2} \
            -S ${rRNA_sam} 2> ${rRNA_log}

    [[ -f ${rRNA_unmap}.1 ]] && \
        mv ${rRNA_unmap}.1 ${rRNA_unmap}_1.fq.gz
    [[ -f ${rRNA_unmap}.2 ]] && \
        mv ${rRNA_unmap}.2 ${rRNA_unmap}_2.fq.gz

    # genome
    local genome_bam="${out_dir}/genome/${fname}_genome.bam"
    local genome_log="${out_dir}/genome/${fname}_genome.log" 
    local genome_unmap="${out_dir}/genome/${fname}_genome_unmap"    
    [[ -f ${bam_q10} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${genome} \
			--un-conc-gz ${genome_unmap} \
            -1 ${rRNA_unmap}_1.fq.gz \
            -2 ${rRNA_unmap}_2.fq.gz 2> ${genome_log} | \
            samtools sort - -O bam -o ${genome_bam} 
        
    [[ -f ${bam_q10} ]] || \
        samtools view ${genome_bam} \
            -q 10 -b -@ 16 -o ${bam_q10}

    [[ -f ${bam_q10}.bai ]] || \
        samtools index ${bam_q10}
    # reporte
    reporte_summary="${out_dir}/${fname}.result"
    if [[ ! -f ${reporte_summary} ]]
    then
	    clean_reads=$(cat ${rRNA_log} | \
			grep "were paired" | awk '{print $1}')
        rRNA_reads=$(cat ${rRNA_log} | \
			grep "1 time" | awk '{sum+=$1}END{print sum}')
        genome_reads=$(cat ${genome_log} | \
			grep "1 time" | awk '{sum+=$1}END{print sum}')
        Unmap_reads=$(cat ${genome_log} | \
			grep "0 time" | awk '{sum+=$1}END{print sum}')
        # add header
        echo "name,clean_reads,rRNA_reads,genome_reads,Unmap_reads" \
			> ${out_dir}/${fname}.stat
        echo "${fname} ${clean_reads} ${rRNA_reads} ${genome_reads} ${Unmap_reads}" | \
			awk '{printf "%s,%d,%d,%d,%d\n",$1,$2,$3,$4,$5}' \
			>> ${out_dir}/${fname}.stat
		# trans to tab
        cat ${out_dir}/${fname}.stat | tr "," "\t" > ${reporte_summary}
		[[ -f ${out_dir}/${fname}.stat ]] && rm ${out_dir}/${fname}.stat
    fi    
    # Remove
    [[ -f ${rRNA_sam} ]] && rm ${rRNA_sam}
    [[ -f ${rRNA_unmap}_1.fq.gz ]] && rm ${rRNA_unmap}_1.fq.gz    
    [[ -f ${rRNA_unmap}_2.fq.gz ]] && rm ${rRNA_unmap}_2.fq.gz
    [[ -f ${genome_bam} ]] && rm ${genome_bam}

}
function cufflinks_FPKM {
	local bam=$1
	local out_dir=$2
	local fname=$(basename ${bam%_genome*})

    FPKM_file="${out_dir}/genes.fpkm_tracking"
    [[ -f ${FPKM_file} ]] || \
        cufflinks -p 10 \
            -o ${out_dir} \
            -G ${gtf_file} \
            ${bam}

    cat ${FPKM_file} | \
        cut -f4,5,10 | tail -n +2 | \
        sort -n -k 2 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\tGene_name\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_FPKM.txt

    Rscript ~/Analysis/qd/src/ID_trans.R \
        ${out_dir} \
        ${fname}_gene_FPKM.txt \
        ${ID_trans_file} \
        ${fname}

    cat ${out_dir}/${fname}_FPKM.txt | tail -n +2 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\tGene_name\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_FPKM.txt

    cat ${out_dir}/${fname}_gene_FPKM.txt | \
        cut -f3 \
        > ${out_dir}/${fname}_gene_FPKM_only.txt

    [[ -f ${out_dir}/${fname}_FPKM.txt ]] && rm ${out_dir}/${fname}_FPKM.txt
}
function bw_covert {
   	local bam=$1
	local out_dir=$2
	local fname=$(basename ${bam%_genome*})

    local bw_file="${out_dir}/${fname}_RPKM.bigWig"
    [[ -f ${bw_file} ]] || \
    bamCoverage -b ${bam} \
        -o ${bw_file} \
        --binSize 30 --smoothLength 90 --numberOfProcessors 15 \
        --extendReads 150 --scaleFactor 1 \
        --normalizeUsing RPKM --centerReads 
} 
function quantify(){
	local bam=$1
    local out_dir=$2
    local gtf=$3
    local fname=$(basename ${bam%_genome*})

	local out_fname="${out_dir}/${fname}.count"
	[[ -f ${out_fname} ]] || \
		featureCounts -T 8 -s 0 \
			-a ${gtf} \
			-o ${out_fname} \
			-M -O --fraction -p -C -B \
			-t exon -g gene_id \
			${bam}

	local count_file="${out_dir}/${fname}_gene_count.txt"
	cat ${out_fname} | tail -n +3 | \
		awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\t"fname1}{print $1"\t"$NF}' \
		> ${count_file}
}
# piepline
fname=$(basename ${fq1%_1.fq*})
# Clean reads
mkdir -p $out_dir/${fname}/01-clean-data

if [ ${method} = "Smart3" ]
then
    Smart3_clean $fq1 $fq2 $out_dir/${fname}/01-clean-data
    wait
elif [ ${method} = "Smart2" ]
then
    Smart2_clean $fq1 $fq2 $out_dir/${fname}/01-clean-data
    wait
else
    echo "No trim adapter"
fi

# Star mapping
mkdir -p $out_dir/${fname}/02-align
align_genome $out_dir/${fname}/01-clean-data/${fname}_1_val_1.fq.gz \
    $out_dir/${fname}/01-clean-data/${fname}_2_val_2.fq.gz \
    $out_dir/${fname}/02-align

# cufflinks counts
mkdir -p $out_dir/${fname}/03-cufflinks
cufflinks_FPKM $out_dir/${fname}/02-align/genome/${fname}_genome_Q10.bam \
    $out_dir/${fname}/03-cufflinks

quantify $out_dir/${fname}/02-align/genome/${fname}_genome_Q10.bam \
	$out_dir/${fname}/03-cufflinks \
	${gtf_file}

mkdir -p $out_dir/${fname}/04-bw
bw_covert $out_dir/${fname}/02-align/genome/${fname}_genome_Q10.bam \
    $out_dir/${fname}/04-bw

# reporte
reporte_summary="${out_dir}/${fname}/${fname}.result"
if [[ ! -f ${reporte_summary} ]]
then
	raw_reads=$(seqkit stat $fq1 | awk '/DNA/{print $4}' | sed 's/,//g')
	clean_reads=$(cat $out_dir/${fname}/02-align/rRNA/${fname}_rRNA.log | \
		grep "were paired" | awk '{print $1}')
	rRNA_reads=$(cat $out_dir/${fname}/02-align/rRNA/${fname}_rRNA.log | \
		grep "1 time" | awk '{sum+=$1}END{print sum}')
	genome_reads=$(cat $out_dir/${fname}/02-align/genome/${fname}_genome.log | \
		grep "1 time" | awk '{sum+=$1}END{print sum}')
	Unmap_reads=$(cat $out_dir/${fname}/02-align/genome/${fname}_genome.log | \
		grep "0 time" | awk '{sum+=$1}END{print sum}')
	# add header
	echo "name,raw_reads,clean_reads,rRNA_reads,genome_reads,Unmap_reads" \
		> ${out_dir}/${fname}.stat
	echo "${fname} ${raw_reads} ${clean_reads} ${rRNA_reads} ${genome_reads} ${Unmap_reads}" | \
		awk '{printf "%s,%d,%d,%d,%d,%d\n",$1,$2,$3,$4,$5,$6}' \
		>> ${out_dir}/${fname}.stat
	# trans to tab
	cat ${out_dir}/${fname}.stat | tr "," "\t" > ${reporte_summary}
	[[ -f ${out_dir}/${fname}.stat ]] && rm ${out_dir}/${fname}.stat
fi    