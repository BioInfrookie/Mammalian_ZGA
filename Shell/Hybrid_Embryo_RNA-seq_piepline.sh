#! /usr/bin/bash

################# Information ###############
# Name : Qiaozhenghao
# Data : 2024-12-23
# Email: qzhrookie@163.com
#############################################

help_str="
参数说明：
    -1) :	fq1
    -2) :	fq2
    -s|--sperm) :	sperm species (Example : mouse; rat; cow; hamster; human; pig)
    -e|--egg) : egg species (Example : mouse; rat; cow; hamster; human; pig)
    -o) :	output directory
    -h) : help document
"

getopt_cmd=$(getopt -o m:e:s:o:1:2:h --long sperm:,egg:,help -n $(basename $0) -- "$@")
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
	    -e|--egg)
            egg_species="$2"
            shift ;;
	    -s|--sperm)
            sperm_species="$2"
            shift ;;
        -m)
            method="$2"
            shift ;;            
        -o)
            out_dir="$2"
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

if [ -n "$sperm_species" ]
then
    if [ ${sperm_species} = "mouse" ]
    then
        sAss="mm10"
        sLabel="mouse"
        sperm_genome="/home/qzh/ref/mm10/hisat2_index/ensemble/genome"
        sperm_rRNA_index="/home/qzh/ref/mm10/hisat2_index/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/mm10/fasta/mm10_ensemble.gtf"
        sperm_ID_trans_file="/home/qzh/ref/mm10/Ensemble_id_trans.txt"
        sTEgtf="/home/qzh/ref/mm10/TEtranscript_ref/mm10_rmsk_TE.gtf"                
    elif [ ${sperm_species} = "cow" ]
    then
        sAss="ARS-UCD1.2"
        sLabel="cow"
        sperm_genome="/home/qzh/ref/ARS-UCD1.2/hisat2_index/genome/genome"
        sperm_rRNA_index="/home/qzh/ref/ARS-UCD1.2/hisat2_index/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/ARS-UCD1.2/fasta/Bos_genome.gtf"
        sperm_ID_trans_file="/home/qzh/ref/ARS-UCD1.2/genome_annotation_data/Ensemble_id_trans.txt"    
        sTEgtf="/home/qzh/ref/ARS-UCD1.2/TEtranscript_ref/bosTau9_rmsk.gtf"
    elif [ ${sperm_species} = "pig" ]
    then
        sAss="ss11"
        sLabel="pig"
        sperm_genome="/home/qzh/ref/ss11/hisat2_index/genome/genome"
        sperm_rRNA_index="/home/qzh/ref/ss11/hisat2_index/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/ss11/fasta/ss11_genome.gtf"
        sperm_ID_trans_file="/home/qzh/ref/ss11/Ensemble_id_trans.txt"
        sTEgtf="/home/qzh/ref/ss11/TEtranscript_ref/susScr11.rmsk.gtf"
    elif [ ${sperm_species} = "Macaca" ]
    then
        sAss="Macaca"
        sLabel="Macaca"
        sperm_genome="/home/qzh/ref/Macaca/hisat2_index/genome/genome"
        sperm_rRNA_index="/home/qzh/ref/Macaca/hisat2_index/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/Macaca/fasta/Macaca_ensemble.gtf"
        sperm_ID_trans_file="/home/qzh/ref/Macaca/annotation/Ensemble_id_trans.txt"
    elif [ ${sperm_species} = "Papio" ]
    then
        sAss="Papio"
        sLabel="Papio"
        sperm_genome="/home/qzh/ref/Papio/hisat2_index/genome/genome"
        sperm_rRNA_index="/home/qzh/ref/Papio/hisat2_index/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/Papio/fasta/Papio_ensemble.gtf"
        sperm_ID_trans_file="/home/qzh/ref/Papio/annotation/Ensemble_id_trans.txt"
    elif [ ${sperm_species} = "human" ]
    then
        sAss="hg38"
        sLabel="human"
		sperm_genome="/home/qzh/ref/hg38/hisat2_index/genome/genome"
		sperm_rRNA_index="/home/qzh/ref/hg38/hisat2_index/rRNA/rRNA"
		sperm_gtf_file="/home/qzh/ref/hg38/fasta/hg38_genome_ensemble.gtf"
		sperm_ID_trans_file="/home/qzh/ref/hg38/genome_annotation_data/Ensemble_id_trans.txt"
        sTEgtf="/home/qzh/ref/hg38/TEtranscript_ref/GRCh38_Ensembl_rmsk_TE.gtf"  
    elif [ ${sperm_species} = "rat" ]
    then
        sAss="rn7"
        sLabel="rat"
        sperm_genome="/home/qzh/ref/rn7/hisat2_index/genome/genome"
        sperm_rRNA_index="/home/qzh/ref/rn7/hisat2_index/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/rn7/genome_annotation_data/rn7_ensemble.gtf"
        sperm_ID_trans_file="/home/qzh/ref/rn7/genome_annotation_data/Ensemble_id_trans.txt"
        sTEgtf="/home/qzh/ref/rn7/TEtranscript_ref/rn7_rmsk_TE.gtf" 
    elif [ ${sperm_species} = "PWK" ]
    then
        sAss="PWK"
        sLabel="PWK"
        sperm_genome="/home/qzh/ref/Mouse_PWK/hisat2/genome/genome"
        sperm_rRNA_index="/home/qzh/ref/Mouse_PWK/hisat2/rRNA/rRNA"
        sperm_gtf_file="/home/qzh/ref/Mouse_PWK/fasta/Mouse_PWK_ensemble.gtf"
        sperm_ID_trans_file="/home/qzh/ref/Mouse_PWK/fasta/Ensemble_id_trans.txt"  
	elif [ ${sperm_species} = "rabbit" ]
	then
	    sAss="rabbit"
        sLabel="rabbit"
		sperm_genome="/home/qzh/ref/OryCun2/hisat2_index/genome/genome"
		sperm_rRNA_index="/home/qzh/ref/OryCun2/hisat2_index/rRNA/rRNA"
		sperm_gtf_file="/home/qzh/ref/OryCun2/fasta/OryCun2_genome_ensemble.gtf"
		sperm_ID_trans_file="/home/qzh/ref/OryCun2/fasta/Ensemble_id_trans.txt"	
	elif [ ${sperm_species} == "hamster" ]
	then
	    sAss="hamster"
        sLabel="hamster"
		sperm_genome="/home/qzh/ref/MesAur1.0/hisat2_index/genome/genome"
		sperm_rRNA_index="/home/qzh/ref/MesAur1.0/hisat2_index/rRNA/rRNA"
		sperm_gtf_file="/home/qzh/ref/MesAur1.0/fasta/MesAur_ensemble_genome.gtf"
		sperm_ID_trans_file="/home/qzh/ref/MesAur1.0/genome_annotation_data/Ensemble_id_trans.txt"	
    else
        echo "species not recognized for sperm"
        exit
    fi
else
    echo "The sperm does not exist."
fi

# prepare egg info
if [ -n "$egg_species" ]
then
    if [ ${egg_species} = "mouse" ]
    then
        eAss="mm10"
        eLabel="mouse"
        egg_genome="/home/qzh/ref/mm10/hisat2_index/ensemble/genome"
        egg_rRNA_index="/home/qzh/ref/mm10/hisat2_index/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/mm10/fasta/mm10_ensemble.gtf"
        egg_ID_trans_file="/home/qzh/ref/mm10/Ensemble_id_trans.txt"    
        eTEgtf="/home/qzh/ref/mm10/TEtranscript_ref/mm10_rmsk_TE.gtf"
    elif [ ${egg_species} = "cow" ]
    then
        eAss="ARS-UCD1.2"
        eLabel="cow"
        egg_genome="/home/qzh/ref/ARS-UCD1.2/hisat2_index/genome/genome"
        egg_rRNA_index="/home/qzh/ref/ARS-UCD1.2/hisat2_index/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/ARS-UCD1.2/fasta/Bos_genome.gtf"
        egg_ID_trans_file="/home/qzh/ref/ARS-UCD1.2/genome_annotation_data/Ensemble_id_trans.txt"    
        eTEgtf="/home/qzh/ref/ARS-UCD1.2/TEtranscript_ref/bosTau9_rmsk.gtf"       
    elif [ ${egg_species} = "pig" ]
    then
        eAss="ss11"
        eLabel="pig"
        egg_genome="/home/qzh/ref/ss11/hisat2_index/genome/genome"
        egg_rRNA_index="/home/qzh/ref/ss11/hisat2_index/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/ss11/fasta/ss11_genome.gtf"
        egg_ID_trans_file="/home/qzh/ref/ss11/Ensemble_id_trans.txt"    
        eTEgtf="/home/qzh/ref/ss11/TEtranscript_ref/susScr11.rmsk.gtf"
    elif [ ${egg_species} = "Macaca" ]
    then
        eAss="Macaca"
        eLabel="Macaca"
        egg_genome="/home/qzh/ref/Macaca/hisat2_index/genome/genome"
        egg_rRNA_index="/home/qzh/ref/Macaca/hisat2_index/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/Macaca/fasta/Macaca_ensemble.gtf"
        egg_ID_trans_file="/home/qzh/ref/Macaca/annotation/Ensemble_id_trans.txt"
    elif [ ${egg_species} = "Papio" ]
    then
        eAss="Papio"
        eLabel="Papio"
        egg_genome="/home/qzh/ref/Papio/hisat2_index/genome/genome"
        egg_rRNA_index="/home/qzh/ref/Papio/hisat2_index/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/Papio/fasta/Papio_ensemble.gtf"
        egg_ID_trans_file="/home/qzh/ref/Papio/annotation/Ensemble_id_trans.txt"   
    elif [ ${egg_species} = "human" ]
    then
        eAss="hg38"
        eLabel="human"
		egg_genome="/home/qzh/ref/hg38/hisat2_index/genome/genome"
		egg_rRNA_index="/home/qzh/ref/hg38/hisat2_index/rRNA/rRNA"
		egg_gtf_file="/home/qzh/ref/hg38/fasta/hg38_genome_ensemble.gtf"
		egg_ID_trans_file="/home/qzh/ref/hg38/genome_annotation_data/Ensemble_id_trans.txt"
        eTEgtf="/home/qzh/ref/hg38/TEtranscript_ref/GRCh38_Ensembl_rmsk_TE.gtf"  
    elif [ ${egg_species} = "rat" ]
    then
        eAss="rn7"
        eLabel="rat"
        egg_genome="/home/qzh/ref/rn7/hisat2_index/genome/genome"
        egg_rRNA_index="/home/qzh/ref/rn7/hisat2_index/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/rn7/genome_annotation_data/rn7_ensemble.gtf"
        egg_ID_trans_file="/home/qzh/ref/rn7/genome_annotation_data/Ensemble_id_trans.txt"
    elif [ ${egg_species} = "PWK" ]
    then
        eAss="PWK"
        eLabel="PWK"
        egg_genome="/home/qzh/ref/Mouse_PWK/hisat2/genome/genome"
        egg_rRNA_index="/home/qzh/ref/Mouse_PWK/hisat2/rRNA/rRNA"
        egg_gtf_file="/home/qzh/ref/Mouse_PWK/fasta/Mouse_PWK_ensemble.gtf"
        egg_ID_trans_file="/home/qzh/ref/Mouse_PWK/fasta/Ensemble_id_trans.txt"
	elif [ ${egg_species} = "rabbit" ]
	then
	    eAss="rabbit"
        eLabel="rabbit"
		egg_genome="/home/qzh/ref/OryCun2/hisat2_index/genome/genome"
		egg_rRNA_index="/home/qzh/ref/OryCun2/hisat2_index/rRNA/rRNA"
		egg_gtf_file="/home/qzh/ref/OryCun2/fasta/OryCun2_genome_ensemble.gtf"
		egg_ID_trans_file="/home/qzh/ref/OryCun2/fasta/Ensemble_id_trans.txt"	
	elif [ ${egg_species} = "hamster" ]
	then
	    eAss="hamster"
        eLabel="hamster"
		egg_genome="/home/qzh/ref/MesAur1.0/hisat2_index/genome/genome"
		egg_rRNA_index="/home/qzh/ref/MesAur1.0/hisat2_index/rRNA/rRNA"
		egg_gtf_file="/home/qzh/ref/MesAur1.0/fasta/MesAur_ensemble_genome.gtf"
		egg_ID_trans_file="/home/qzh/ref/MesAur1.0/genome_annotation_data/Ensemble_id_trans.txt"		
    else
        echo "species not recognized for egg"
        exit
    fi
else
    echo "The egg does not exist."
fi

echo "####### perform analysis ######"
echo "species of sperm: $sLabel"
echo "assembly of sperm: $sAss"
echo "species of egg: $eLabel"
echo "assembly of egg: $eAss"

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
        cutadapt -j 20 -n 2 -m 40 -O 5 -e 0.2 \
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
			-a A{10} -a T{10} \
			-a2 A{10} -a2 T{10} \
            --paired --clip_R1 10 --clip_R2 10 \
            --three_prime_clip_R1 3 --three_prime_clip_R2 3 \
            ${fq1} ${fq2} \
            -o ${out_dir}
} 
function Select_fq() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1_*})

    local fq1_both_map="${out_dir}/${fname}_both_1.fq.gz"
    local fq2_both_map="${out_dir}/${fname}_both_2.fq.gz" 

    # rRNA
    local egg_rRNA_sam="${out_dir}/${fname}_egg_rRNA.sam"
    local egg_rRNA_unmap="${out_dir}/${fname}_egg_rRNA_unmap"
    [[ -f ${fq1_both_map} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${egg_rRNA_index} \
            --un-conc-gz ${egg_rRNA_unmap} \
            -1 ${fq1} \
            -2 ${fq2} \
            -S ${egg_rRNA_sam}

    [[ -f ${egg_rRNA_unmap}.1 ]] && \
        mv ${egg_rRNA_unmap}.1 ${egg_rRNA_unmap}_1.fq.gz
    [[ -f ${egg_rRNA_unmap}.2 ]] && \
        mv ${egg_rRNA_unmap}.2 ${egg_rRNA_unmap}_2.fq.gz

    local sperm_rRNA_sam="${out_dir}/${fname}_sperm_rRNA.sam"
    local sperm_rRNA_unmap="${out_dir}/${fname}_sperm_rRNA_unmap"
    [[ -f ${fq1_both_map} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${sperm_rRNA_index} \
            --un-conc-gz ${sperm_rRNA_unmap} \
            -1 ${egg_rRNA_unmap}_1.fq.gz \
            -2 ${egg_rRNA_unmap}_2.fq.gz \
            -S ${sperm_rRNA_sam}

    [[ -f ${sperm_rRNA_unmap}.1 ]] && \
        mv ${sperm_rRNA_unmap}.1 ${sperm_rRNA_unmap}_1.fq.gz
    [[ -f ${sperm_rRNA_unmap}.2 ]] && \
        mv ${sperm_rRNA_unmap}.2 ${sperm_rRNA_unmap}_2.fq.gz

    # map to sperm
    local sperm_genome_sam="${out_dir}/${fname}_to_${sLabel}.sam"
    local sperm_genome_unmap="${out_dir}/${fname}_to_${sLabel}_unmap"
    [[ -f ${fq1_both_map} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${sperm_genome} \
            --un-conc-gz ${sperm_genome_unmap} \
            -1 ${sperm_rRNA_unmap}_1.fq.gz \
            -2 ${sperm_rRNA_unmap}_2.fq.gz \
            -S ${sperm_genome_sam} 
    [[ -f ${sperm_genome_unmap}.1 ]] && \
        mv ${sperm_genome_unmap}.1 ${sperm_genome_unmap}_1.fq.gz
    [[ -f ${sperm_genome_unmap}.2 ]] && \
        mv ${sperm_genome_unmap}.2 ${sperm_genome_unmap}_2.fq.gz

    # map to egg
    local egg_genome_sam="${out_dir}/${fname}_to_${eLabel}.sam"
	local both_unmap="${out_dir}/${fname}_both_unmap"
    [[ -f ${fq1_both_map} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
			--un-conc-gz ${both_unmap} \
            -x ${egg_genome} \
            -1 ${sperm_genome_unmap}_1.fq.gz \
            -2 ${sperm_genome_unmap}_2.fq.gz \
            -S ${egg_genome_sam} 

    [[ -f ${both_unmap}.1 ]] && \
        mv ${both_unmap}.1 ${both_unmap}_1.fq.gz
    [[ -f ${both_unmap}.2 ]] && \
        mv ${both_unmap}.2 ${both_unmap}_2.fq.gz

    # with-sperm map to egg
    local fq1_map_sperm="${out_dir}/${fname}_to_${sLabel}_1.fq.gz"
    local fq2_map_sperm="${out_dir}/${fname}_to_${sLabel}_2.fq.gz"
    [[ -f ${fq1_both_map} ]] || \
        samtools fastq -F 12 -@ 10 \
            -1 ${fq1_map_sperm} \
            -2 ${fq2_map_sperm} \
            ${sperm_genome_sam} 

    local both_genome_sam="${out_dir}/${fname}_both.sam"
    local only_sperm="${out_dir}/${fname}_only_to_${sLabel}"
    [[ -f ${fq1_both_map} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${egg_genome} \
            --un-conc-gz ${only_sperm} \
            -1 ${fq1_map_sperm} \
            -2 ${fq2_map_sperm} \
            -S ${both_genome_sam} 

    [[ -f ${only_sperm}.1 ]] && \
        mv ${only_sperm}.1 ${only_sperm}_1.fq.gz
    [[ -f ${only_sperm}.2 ]] && \
        mv ${only_sperm}.2 ${only_sperm}_2.fq.gz

    # only map to egg
    local fq1_only_egg="${out_dir}/${fname}_only_to_${eLabel}_1.fq.gz"
    local fq2_only_egg="${out_dir}/${fname}_only_to_${eLabel}_2.fq.gz"
    [[ -f ${fq1_only_egg} ]] || \
        samtools fastq -F 12 -@ 10 \
            -1 ${fq1_only_egg} \
            -2 ${fq2_only_egg} \
            ${egg_genome_sam}

    # both map to sperm and egg
    [[ -f ${fq1_both_map} ]] || \
        samtools fastq -F 12 -@ 10 \
            -1 ${fq1_both_map} \
            -2 ${fq2_both_map} \
            ${both_genome_sam}

    # Remove
    [[ -f ${fq1_map_sperm} ]] && rm ${fq1_map_sperm}
    [[ -f ${fq2_map_sperm} ]] && rm ${fq2_map_sperm}
    [[ -f ${both_genome_sam} ]] && rm ${out_dir}/*sam
    [[ -f ${egg_rRNA_unmap}_1.fq.gz ]] && rm ${egg_rRNA_unmap}_1.fq.gz ${egg_rRNA_unmap}_2.fq.gz 
    # [[ -f ${sperm_rRNA_unmap}_1.fq.gz ]] && rm ${sperm_rRNA_unmap}_1.fq.gz ${sperm_rRNA_unmap}_2.fq.gz 
    [[ -f ${sperm_genome_unmap}_1.fq.gz ]] && rm ${sperm_genome_unmap}_1.fq.gz ${sperm_genome_unmap}_2.fq.gz
	[[ -f ${Human_genome_bam}  ]] && rm ${Human_genome_bam} 
	
	
}
function Align_STAR_egg() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

    local genome_bam="${out_dir}/${fname}.bam" 
    local genome_bam_q10="${out_dir}/${fname}_Q10.bam"      
    [[ -f ${genome_bam_q10} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${egg_genome} \
            -1 ${fq1} \
            -2 ${fq2} | \
            samtools sort - -O bam -o ${genome_bam} 
        
    [[ -f ${genome_bam_q10} ]] || \
        samtools view ${genome_bam} \
            -q 10 -b -@ 16 -o ${genome_bam_q10}

    [[ -f ${genome_bam_q10}.bai ]] || \
        samtools index ${genome_bam_q10}  
    # Remove
    [[ -f ${genome_bam} ]] && rm ${genome_bam}

}
function Align_STAR_sperm() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

    local genome_bam="${out_dir}/${fname}.bam"
    local genome_bam_q10="${out_dir}/${fname}_Q10.bam"      
    [[ -f ${genome_bam_q10} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${sperm_genome} \
            -1 ${fq1} \
            -2 ${fq2} | \
            samtools sort - -O bam -o ${genome_bam} 
        
    [[ -f ${genome_bam_q10} ]] || \
        samtools view ${genome_bam} \
            -q 10 -b -@ 16 -o ${genome_bam_q10}

    [[ -f ${genome_bam_q10}.bai ]] || \
        samtools index ${genome_bam_q10}
    # Remove
    [[ -f ${genome_bam} ]] && rm ${genome_bam}
}
function Align_STAR_both() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

    local sperm_align_bam="${out_dir}/${fname}_to_${sLabel}.bam"
    local sperm_align_q10="${out_dir}/${fname}_to_${sLabel}_Q10.bam"      
    [[ -f ${sperm_align_q10} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${sperm_genome} \
            -1 ${fq1} \
            -2 ${fq2} | \
            samtools sort - -O bam -o ${sperm_align_bam} 
        
    [[ -f ${sperm_align_q10} ]] || \
        samtools view ${sperm_align_bam} \
            -q 10 -b -@ 16 -o ${sperm_align_q10}

    local egg_align_bam="${out_dir}/${fname}_to_${eLabel}.bam" 
    local egg_align_q10="${out_dir}/${fname}_to_${eLabel}_Q10.bam" 
    [[ -f ${egg_align_q10} ]] || \
        hisat2 -p 15 --dta-cufflinks \
            --no-mixed --no-discordant \
            -x ${egg_genome} \
            -1 ${fq1} \
            -2 ${fq2} | \
            samtools sort - -O bam -o ${egg_align_bam} 
        
    [[ -f ${egg_align_q10} ]] || \
        samtools view ${egg_align_bam} \
            -q 10 -b -@ 16 -o ${egg_align_q10}

    # Remove
    [[ -f ${sperm_align_bam} ]] && rm ${sperm_align_bam}
    [[ -f ${egg_align_bam} ]] && rm ${egg_align_bam}   
}
function Gene_Count_egg() {
    local bam=$1
    local out_dir=$2
    local fname=$(basename ${bam%_Q10*})

    FPKM_file="${out_dir}/genes.fpkm_tracking"
    [[ -f ${FPKM_file} ]] || \
        cufflinks -p 10 \
            -o ${out_dir} \
            -G ${egg_gtf_file} \
            ${bam}

    cat ${FPKM_file} | \
        cut -f4,5,10 | tail -n +2 | \
        sort -n -k 2 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\tGene_name\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_FPKM.txt

    Rscript ~/Analysis/qd/src/ID_trans.R \
        ${out_dir} \
        ${fname}_gene_FPKM.txt \
        ${egg_ID_trans_file} \
        ${fname}

    cat ${out_dir}/${fname}_FPKM.txt | tail -n +2 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\tGene_name\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_FPKM.txt

    cat ${out_dir}/${fname}_gene_FPKM.txt | \
        cut -f3 \
        > ${out_dir}/${fname}_gene_FPKM_only.txt

    [[ -f ${out_dir}/${fname}_FPKM.txt ]] && rm ${out_dir}/${fname}_FPKM.txt

}
function Gene_Count_sperm() {
    local bam=$1
    local out_dir=$2
    local fname=$(basename ${bam%_Q10*})

    FPKM_file="${out_dir}/genes.fpkm_tracking"
    [[ -f ${FPKM_file} ]] || \
        cufflinks -p 10 \
            -o ${out_dir} \
            -G ${sperm_gtf_file} \
            ${bam}

    cat ${FPKM_file} | \
        cut -f4,5,10 | tail -n +2 | \
        sort -n -k 2 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\tGene_name\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_FPKM.txt

    Rscript ~/Analysis/qd/src/ID_trans.R \
        ${out_dir} \
        ${fname}_gene_FPKM.txt \
        ${sperm_ID_trans_file} \
        ${fname}

    cat ${out_dir}/${fname}_FPKM.txt | tail -n +2 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\tGene_name\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_FPKM.txt

    cat ${out_dir}/${fname}_gene_FPKM.txt | \
        cut -f3 \
        > ${out_dir}/${fname}_gene_FPKM_only.txt

    [[ -f ${out_dir}/${fname}_FPKM.txt ]] && rm ${out_dir}/${fname}_FPKM.txt
}
function TE_RPM_egg() {
    local bam=$1
    local out_dir=$2
    local fname=$(basename ${bam%_Aligned*})

    local egg_TE="${out_dir}/${fname}_TE_count.cntTable"
    [[ -f ${egg_TE} ]] || \
        TEcount -b ${bam} \
                --GTF ${egtf} \
                --TE ${eTEgtf} \
                --mode multi \
                --sortByPos \
                --project ${fname}_TE_count \
                --outdir ${out_dir}

    local egg_TE_count="${out_dir}/${fname}_TE_count.txt"
    cat ${egg_TE} | tail -n +2 | grep -v "^\"EN" | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "id\t"fname1}{print $0}' \
        > ${egg_TE_count}

    local egg_TE_RPM="${out_dir}/${fname}_TE_RPM.txt"
    cat ${egg_TE_count} | \
        awk '{print $1"\t"$2*"'"${scale_facor}"'"}' \
        > ${egg_TE_RPM}

    local egg_TE_RPM_only="${out_dir}/${fname}_TE_RPM_only.txt"
    cat ${egg_TE_RPM} | cut -f2 | awk -v fname1="'"${fname}"'" 'BEGIN{print fname1}{print $0}' \
        > ${egg_TE_RPM_only}
}
function TE_RPM_sperm() {
    local bam=$1
    local out_dir=$2
    local fname=$(basename ${bam%_Aligned*})

    local sperm_TE="${out_dir}/${fname}_TE_count.cntTable"
    [[ -f ${sperm_TE} ]] || \
        TEcount -b ${bam} \
                --GTF ${sgtf} \
                --TE ${sTEgtf} \
                --mode multi \
                --sortByPos \
                --project ${fname}_TE_count \
                --outdir ${out_dir}

    local sperm_TE_count="${out_dir}/${fname}_TE_count.txt"
    cat ${sperm_TE} | tail -n +2 | grep -v "^\"EN" | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "id\t"fname1}{print $0}' \
        > ${sperm_TE_count}

    local sperm_TE_RPM="${out_dir}/${fname}_TE_RPM.txt"
    cat ${sperm_TE_count} | \
        awk '{print $1"\t"$2*"'"${scale_facor}"'"}' \
        > ${sperm_TE_RPM}

    local sperm_TE_RPM_only="${out_dir}/${fname}_TE_RPM_only.txt"
    cat ${sperm_TE_RPM} | cut -f2 | awk -v fname1="'"${fname}"'" 'BEGIN{print fname1}{print $0}' \
        > ${sperm_TE_RPM_only}
}
function bw_covert {
   	local bam=$1
	local out_dir=$2
	local fname=$(basename ${bam%_Q10*})

    local bw_file="${out_dir}/${fname}_RPKM.bigWig"
    [[ -f ${bw_file} ]] || \
		bamCoverage -b ${bam} \
			-o ${bw_file} \
			--binSize 30 --smoothLength 90 --numberOfProcessors 15 \
			--extendReads 150 --scaleFactor 1 \
			--normalizeUsing RPKM --centerReads    
} 

# piepline 
fname=$(basename ${fq1%_1.fq*})
mkdir -p ${out_dir}/${fname}
## 01-clean
mkdir -p ${out_dir}/${fname}/01-clean
if [ ${method} = "Smart3" ]
then
    Smart3_clean $fq1 $fq2 ${out_dir}/${fname}/01-clean
    wait
elif [ ${method} = "Smart2" ]
then
    Smart2_clean $fq1 $fq2 ${out_dir}/${fname}/01-clean
    wait
else
    echo "No trim adapter"
fi

## 02-select
mkdir -p ${out_dir}/${fname}/02-Select_egg_sperm_fq
Select_fq ${out_dir}/${fname}/01-clean/${fname}_1_val_1.fq.gz \
    ${out_dir}/${fname}/01-clean/${fname}_2_val_2.fq.gz \
    ${out_dir}/${fname}/02-Select_egg_sperm_fq

## 03-align
mkdir -p ${out_dir}/${fname}/03-STAR_align
Align_STAR_egg ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_only_to_${eLabel}_1.fq.gz \
    ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_only_to_${eLabel}_2.fq.gz \
    ${out_dir}/${fname}/03-STAR_align

Align_STAR_sperm ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_only_to_${sLabel}_1.fq.gz \
    ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_only_to_${sLabel}_2.fq.gz \
    ${out_dir}/${fname}/03-STAR_align

Align_STAR_both ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_both_1.fq.gz \
    ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_both_2.fq.gz \
    ${out_dir}/${fname}/03-STAR_align

## 04-gene-count
mkdir -p ${out_dir}/${fname}/04-gene_count/only_sperm
Gene_Count_sperm ${out_dir}/${fname}/03-STAR_align/${fname}_only_to_${sLabel}_Q10.bam \
    ${out_dir}/${fname}/04-gene_count/only_sperm

mkdir -p ${out_dir}/${fname}/04-gene_count/only_egg
Gene_Count_egg ${out_dir}/${fname}/03-STAR_align/${fname}_only_to_${eLabel}_Q10.bam \
    ${out_dir}/${fname}/04-gene_count/only_egg

mkdir -p ${out_dir}/${fname}/04-gene_count/both_egg_sperm
Gene_Count_egg ${out_dir}/${fname}/03-STAR_align/${fname}_both_to_${eLabel}_Q10.bam \
    ${out_dir}/${fname}/04-gene_count/both_egg_sperm

Gene_Count_sperm ${out_dir}/${fname}/03-STAR_align/${fname}_both_to_${sLabel}_Q10.bam \
    ${out_dir}/${fname}/04-gene_count/both_egg_sperm

## 05-bw
mkdir -p $out_dir/${fname}/05-bw
mkdir -p ${out_dir}/${fname}/05-bw/only_sperm
bw_covert ${out_dir}/${fname}/03-STAR_align/${fname}_only_to_${sLabel}_Q10.bam \
    ${out_dir}/${fname}/05-bw/only_sperm

mkdir -p ${out_dir}/${fname}/05-bw/only_egg
bw_covert ${out_dir}/${fname}/03-STAR_align/${fname}_only_to_${eLabel}_Q10.bam \
    ${out_dir}/${fname}/05-bw/only_egg


report_file="${out_dir}/${fname}/${fname}_reads.result"
if [[ ! -f ${report_file} ]]
then
    Total_reads=`seqkit stat ${fq1} | \
        awk '/DNA/{print $4}' | sed 's/,//g'`
    Clean_reads=`seqkit stat ${out_dir}/${fname}/01-clean/${fname}_1_val_1.fq.gz | \
        awk '/DNA/{print $4}' | sed 's/,//g'`
    Only_sperm_reads=`seqkit stat ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_only_to_${sLabel}_1.fq.gz | \
        awk '/DNA/{print $4}' | sed 's/,//g'`
    Only_egg_reads=`seqkit stat ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_only_to_${eLabel}_1.fq.gz | \
        awk '/DNA/{print $4}' | sed 's/,//g'`
    Both_sperm_egg_reads=`seqkit stat ${out_dir}/${fname}/02-Select_egg_sperm_fq/${fname}_both_1.fq.gz | \
        awk '/DNA/{print $4}' | sed 's/,//g'`
    Unmap_reads=`echo ${Clean_reads} ${Only_sperm_reads} ${Only_egg_reads} ${Both_sperm_egg_reads} | \
        awk '{print $1-$2-$3-$4}'`
    Only_sperm_gene_num=`cat ${out_dir}/${fname}/04-gene_count/only_sperm/${fname}_only_to_${sLabel}_gene_FPKM_only.txt | \
        awk '{if($1>0)print $0}'| wc -l `
    Only_egg_gene_num=`cat ${out_dir}/${fname}/04-gene_count/only_egg/${fname}_only_to_${eLabel}_gene_FPKM_only.txt | \
        awk '{if($1>0)print $0}'| wc -l `
    
    echo "Sample,Total_reads,Clean_reads,Both_sperm_egg_reads,Only_egg_reads,Only_sperm_reads,Unmap_reads,Only_sperm_gene_num,Only_egg_gene_num" \
        > ${out_dir}/${fname}/${fname}.stat
    echo "${fname} ${Total_reads} ${Clean_reads} ${Both_sperm_egg_reads} ${Only_egg_reads} ${Only_sperm_reads} ${Unmap_reads} ${Only_sperm_gene_num} ${Only_egg_gene_num}" | \
        awk '{printf "%s,%d,%d,%d,%d,%d,%d,%d,%d\n",$1,$2,$3,$4,$5,$6,$7,$8,$9}' \
        >> ${out_dir}/${fname}/${fname}.stat
    # trans to tab
    cat ${out_dir}/${fname}/${fname}.stat | tr "," "\t" > ${report_file}

	[[ -f ${out_dir}/${fname}/${fname}.stat ]] && rm ${out_dir}/${fname}/${fname}.stat
fi

