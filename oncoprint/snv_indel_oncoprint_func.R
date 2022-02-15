library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

process_metadata = function(purity_file="absolute_purity_ploidy.txt", patient_file="patients.txt"){
  # load purity and clinical information
  purity = read.table(purity_file, header=TRUE, as.is=TRUE)
  patients = read.table(patient_file, header=TRUE, as.is=TRUE)[,c("patient","clin_ben","pfs_mon")]

  # create metadata by merging purity and patient information
  meta = merge(purity, patients, by=c("patient"))
  rownames(meta) = meta$sample_name
  meta$tumor_fraction = rep("Low", nrow(meta))
  meta[meta$absolute_purity >= 0.15, "tumor_fraction"] = "High"
  
  return(meta)
}

prepare_oncopanel_scheme = function(){
  col = c("springgreen4", "darkorchid1", "blue4", "yellow1", "darkorange", "red4")
  names(col) = c("missense", "splice_site", "nonsense", "frame_shift", "in_frame_indel", "other_non_syn")

  alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # other non-synonymous mutation
    other_non_syn = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["other_non_syn"], col = NA))
    },
    # frame shift
    frame_shift = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["frame_shift"], col = NA))
    },
    # in frame indel
    in_frame_indel = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["in_frame_indel"], col = NA))
    },
    # missense
    missense = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["missense"], col = NA))
    },
    # nonsense
    nonsense = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["nonsense"], col = NA))
    },
    # splice site
    splice_site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
            gp = gpar(fill = col["splice_site"], col = NA))
    }
  )
  return(list(alter_fun=alter_fun, col=col))
}

load_matrix = function(matrix_file, metadata, gene_list=c(NA)){
  mat = read.table(matrix_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)
  mat[is.na(mat)] = ""
  if(!is.na(gene_list[1])){
    mat = mat[gene_list, ]
  }
  meta = metadata[colnames(mat),]
  return(list(mat=mat, meta=meta))
}

plot_onco_print = function(data, onco_setting, legend_param, row_order, column_order, title="", remove_empty_rows=TRUE, patient=NA){
  gradient_fun = colorRamp2(c(2, 4), c("white", "purple"))
  gradient_fun2 = colorRamp2(c(1, 32), c("black", "red"))
  annot_col = list( tm=c("Baseline"="firebrick3", "Cycle1"="gold", "EOT"="dodgerblue4", "Tissue"="springgreen"),
                  tf=c("High"="darkslategray3", "Low"="darkslategray1"),
                  pl=gradient_fun,
                  cb=c("Yes"="steelblue", "No"="tomato"))
  if(!is.na(patient)){
      ppl = brewer.pal("Set3",n=length(patient))
      names(ppl) = as.character(patient)
      annot_col = list( tm=c("Baseline"="firebrick3", "Cycle1"="gold", "EOT"="dodgerblue4", "Tissue"="springgreen"),
                      tf=c("High"="darkslategray3", "Low"="darkslategray1"),
                      pl=gradient_fun,
                      cb=c("Yes"="steelblue", "No"="tomato"),
                      pp=ppl)
  }
  onco_t = oncoPrint(data$mat, alter_fun = onco_setting$alter_fun, col = onco_setting$col, column_title=title,
                     heatmap_legend_param=legend_param,
                     remove_empty_rows = remove_empty_rows,
                     pct_side = "right", row_names_side = "left",
                     right_annotation = rowAnnotation(
                     row_barplot = anno_oncoprint_barplot(
                         type=c("missense", "splice_site", "nonsense", "frame_shift", "in_frame_indel", "other_non_syn"))),
                     top_annotation = HeatmapAnnotation(
                     tm=data$meta[colnames(data$mat),"timepoint"],
                     tf=data$meta[colnames(data$mat),"tumor_fraction"],
                     pl=data$meta[colnames(data$mat),"absolute_ploidy"],
                     cb=data$meta[colnames(data$mat),"clin_ben"],
                     pp=data$meta[colnames(data$mat),"patient"],
                     col=annot_col))
  if(is.na(column_order)){
    clin_ben_strs = data$meta$clin_ben
    column_order = c(onco_t@column_order[which(clin_ben_strs[onco_t@column_order] == "No")],
                     onco_t@column_order[which(clin_ben_strs[onco_t@column_order] == "Yes")])
  }
  return( oncoPrint(data$mat, alter_fun = onco_setting$alter_fun, col = onco_setting$col, column_title=title,
                    heatmap_legend_param=legend_param,
                    column_order=column_order,
                    row_order=row_order,
                    show_column_names=TRUE,
                    remove_empty_rows = remove_empty_rows,
                    pct_side = "right", row_names_side = "left",
                    pct_gp=gpar(fontsize=6),
                    column_names_gp=gpar(fontsize=5),
                    right_annotation = rowAnnotation(
                    row_barplot = anno_oncoprint_barplot(
                           type=c("missense", "splice_site", "nonsense", "frame_shift", "in_frame_indel", "other_non_syn"))),
                    top_annotation = HeatmapAnnotation(
                    tm=data$meta[colnames(data$mat),"timepoint"],
                    tf=data$meta[colnames(data$mat),"tumor_fraction"],
                    pl=data$meta[colnames(data$mat),"absolute_ploidy"],
                    cb=data$meta[colnames(data$mat),"clin_ben"],
                    pp=data$meta[colnames(data$mat),"patient"],
                    col=annot_col, simple_anno_size = unit(2, "mm"))) )
}

