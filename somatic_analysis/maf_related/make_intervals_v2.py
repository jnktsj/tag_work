import sys,collections

sample_id=sys.argv[1]
#this_maf=sys.argv[2]
print "reading:" +sample_id
all_variants={}
cols_to_drop=set(['t_alt_count','t_ref_count','i_ffpe_p_value', 'i_oxog_p_value', 'i_t_ALT_F1R2', 'i_t_ALT_F2R1', 'i_t_Foxog', 'i_t_REF_F1R2', 'i_t_REF_F2R1', 'id', 'init_n_lod', 'init_t_lod', 'map_q0_reads', 'n_lod', 'n_q20_count', 'normal_power', 'normal_power_nsp', 'normal_power_wsp', 'phasing_genotype', 'phasing_id', 'phred_scaled_likelihoods', 'power', 'power_to_detect_negative_strand_artifact', 'power_to_detect_positive_strand_artifact', 'qual', 'read_depth', 'refseq_mrna_id', 'repeat_times_tandem_repeat_unit', 'score', 'secondary_variant_classification', 'set', 'short_tandem_repeat_membership', 'strand_bias_counts', 't_alt_max_mapq', 't_del_count', 't_ins_count', 't_lod_fstar', 't_lod_fstar_forward', 't_lod_fstar_reverse', 't_q20_count', 't_ref_max_mapq', 'total_reads', 'tumor_alt_fpir_mad', 'tumor_alt_fpir_median', 'tumor_alt_rpir_mad', 'tumor_alt_rpir_median', 'tumor_f', 'tumor_power'])

for tsv_fn in sys.argv[2:]:
    print "loading:" +tsv_fn
    file_handle=open(tsv_fn,'rU')
    header=file_handle.readline()
    while header[0] == "#" or not header.strip():
        header=file_handle.readline()

    header=header.strip("\n").split("\t")
    h=collections.OrderedDict([[x[1],x[0]] for x in enumerate(header)])
    for line in file_handle:
        fields=line.strip("\n").split("\t")
        if "filter" in h and fields[h["filter"]] != "PASS":continue
        var_str=":".join([fields[h["Chromosome"]],fields[h["Start_position"]],fields[h["Tumor_Seq_Allele2"]]])
        if var_str not in all_variants:
            all_variants[var_str]=collections.OrderedDict([[x[1],fields[x[0]]] for x in enumerate(header) if x[1] not in cols_to_drop])
    file_handle.close()

#write interval file
maf_variants={}
print "writing:" + sample_id

for maf_fn in sys.argv[2:]:
    if sample_id in maf_fn:this_maf=maf_fn #Infer the correct maf from filenames, because firecloud cannot handle the sample input twice 

file_handle=open(this_maf,'rU')
header=file_handle.readline()
while header[0] == "#" or not header.strip():
    header=file_handle.readline()

header=header.strip("\n").split("\t")
h=collections.OrderedDict([[x[1],x[0]] for x in enumerate(header)])
for line in file_handle:
    fields=line.strip("\n").split("\t")
    var_str=":".join([fields[h["Chromosome"]],fields[h["Start_position"]],fields[h["Tumor_Seq_Allele2"]]])
    maf_variants[var_str]=collections.OrderedDict([[x[1],fields[x[0]]] for x in enumerate(header)])

with open(sample_id+".forecallready.maf","w") as out_maf:
    out_maf.write("\t".join(maf_variants.values()[0].keys())+"\n")
    for variant in all_variants:
        if variant in maf_variants:
            out_maf.write("\t".join([maf_variants[variant][x] for x in header])+"\n")
        else:
            out_maf.write("\t".join([all_variants[variant][x] if x in all_variants[variant] else "" for x in header])+"\n")

