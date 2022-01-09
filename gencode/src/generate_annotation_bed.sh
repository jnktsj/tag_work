# coding regions restricted with CCDS and basic tags
# !! many known transcripts do not have the CCDS labels, so don't use this !!
# grep 'tag "CCDS"' ../gencode.v19.annotation.gtf | grep 'tag "basic"' | \
# awk 'BEGIN{OFS="\t"; FS="\t"}
#     {split($1,chrom,"chr"); split($9,raw,"gene_name "); split(raw[2],gene,";"); split($9,raw,"gene_type "); split(raw[2],type,";");
#      if($3 == "CDS"){if(chrom[2] == "M"){chrom[2]="MT"} print chrom[2],$4-1,$5,gene[1]":"type[1]":"$7}}' | \
# sed 's/"//g' | sort -k4,4 -k1,1 -k2,2n | uniq | python ../src/merge_entries.py - > gencode_v19_CDS.CCDS.bed

# CDS restricted with basic tag
grep 'tag "basic"' ../gencode.v19.annotation.gtf | \
awk 'BEGIN{OFS="\t"; FS="\t"}
     {split($1,chrom,"chr"); split($9,raw,"gene_name "); split(raw[2],gene,";"); split($9,raw,"gene_type "); split(raw[2],type,";");
      if($3 == "CDS"){if(chrom[2] == "M"){chrom[2]="MT"} print chrom[2],$4-1,$5,gene[1]":"type[1]":"$7}}' | \
sed 's/"//g' | sort -k4,4 -k1,1 -k2,2n | uniq | python ../src/merge_entries.py - > gencode_v19_CDS.bed

# Exon restricted with basic tag
grep 'tag "basic"' ../gencode.v19.annotation.gtf | \
awk 'BEGIN{OFS="\t"; FS="\t"}
     {split($1,chrom,"chr"); split($9,raw,"gene_name "); split(raw[2],gene,";"); split($9,raw,"gene_type "); split(raw[2],type,";");
      if($3 == "exon"){if(chrom[2] == "M"){chrom[2]="MT"} print chrom[2],$4-1,$5,gene[1]":"type[1]":"$7}}' | \
sed 's/"//g' | sort -k4,4 -k1,1 -k2,2n | uniq | python ../src/merge_entries.py - > gencode_v19_exon.bed

# Intron, TSS/TES
python ../src/get_intron_TSS_TES.py ../hg19_length.txt gencode_v19_exon.bed 
mv TSS_TES.bed gencode_v19_TSS-TES.bed
mv introns.bed gencode_v19_intron.bed

# Format CDS and Exon
awk '{print $0":CDS"}' gencode_v19_CDS.bed > tmp
mv tmp gencode_v19_CDS.bed
awk '{print $0":Exon"}' gencode_v19_exon.bed > tmp
mv tmp gencode_v19_exon.bed


# Create non-overlapping gene model
bedtools subtract -a gencode_v19_exon.bed -b gencode_v19_CDS.bed > non_coding_exon.bed
bedtools subtract -a gencode_v19_TSS-TES.bed -b gencode_v19_exon.bed > pure_TSS-TES.bed
cat gencode_v19_exon.bed gencode_v19_TSS-TES.bed | bedtools sort -i stdin | bedtools merge -i stdin | bedtools subtract -a gencode_v19_intron.bed -b stdin > pure_intron.bed
cat gencode_v19_*bed | bedtools sort -i stdin | bedtools merge -i stdin | bedtools subtract -a ../hg19.bed -b stdin | awk '{print $0"\tIntergenic"}' > intergenic.bed
cat gencode_v19_CDS.bed non_coding_exon.bed pure_TSS-TES.bed pure_intron.bed intergenic.bed | bedtools sort -i stdin > hg19_annotation.gencode_v19_gene_model.bed
rm intergenic.bed non_coding_exon.bed pure_TSS-TES.bed pure_intron.bed

=====


# print CDS overlap with ICE interval list
echo "========="
echo "CDS overlap in ICE baits"
grep -v ^@ /seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list | \
awk '{print $1,$2-1,$3}' OFS="\t" | bedtools intersect -a gene_list_cds_merged.bed -b stdin | \
sort -k1,1 -k2,2n -k4,4 | python merge_entries.py - | \
awk '{arr[$4] += $3-$2}END{for(i in arr){print i"\t"arr[i]}}' | \
paste <(awk '{arr[$4] += $3-$2}END{for(i in arr){print i"\t"arr[i]}}' gene_list_cds_merged.bed) - | awk '{print $1,$2,$4,$4/$2}' OFS="\t"
echo "---------"
echo "CDS overlap in ICE targets"
grep -v ^@ /seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list | \
awk '{print $1,$2-1,$3}' OFS="\t" | bedtools intersect -a gene_list_cds_merged.bed -b stdin | \
sort -k1,1 -k2,2n -k4,4 | python merge_entries.py - | \
awk '{arr[$4] += $3-$2}END{for(i in arr){print i"\t"arr[i]}}' | \
paste <(awk '{arr[$4] += $3-$2}END{for(i in arr){print i"\t"arr[i]}}' gene_list_cds_merged.bed) - | awk '{print $1,$2,$4,$4/$2}' OFS="\t"



# print CDS overlap with Twist interval list
echo "========="
echo "CDS overlap in Twist probes"
grep -v ^@ /broad/tag_working/jtsuji/tag_263/Twist/Exome_V_1_3_0_Probes_hg19.liftedOver.interval_list | \
awk '{print $1,$2-1,$3}' OFS="\t" | bedtools intersect -a gene_list_cds_merged.bed -b stdin | \
sort -k1,1 -k2,2n -k4,4 | python merge_entries.py - | \
awk '{arr[$4] += $3-$2}END{for(i in arr){print i"\t"arr[i]}}' | \
paste <(awk '{arr[$4] += $3-$2}END{for(i in arr){print i"\t"arr[i]}}' gene_list_cds_merged.bed) - | awk '{print $1,$2,$4,$4/$2}' OFS="\t"
