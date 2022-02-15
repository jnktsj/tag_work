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
  col = c("red4", "hotpink1", "steelblue2", "blue3", "chartreuse3")
  names(col) = c("amplification", "gain", "loss", "deletion", "nLOH")

  alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # copy number LOH
    nLOH = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["nLOH"], col = NA))
    },
    # deletion
    deletion = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["deletion"], col = NA))
    },
    # loss
    loss = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["loss"], col = NA))
    },
    # gain
    gain = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["gain"], col = NA))
    },
    # amplification
    amplification = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
            gp = gpar(fill = col["amplification"], col = NA))
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

plot_onco_print = function(data, onco_setting, legend_param, row_order=c(NA), column_order=c(NA), title="", remove_empty_rows=TRUE){
  gradient_fun = colorRamp2(c(2, 4), c("white", "purple"))
  gradient_fun2 = colorRamp2(c(1, 32), c("black", "red"))
  annot_col = list( tm=c("Baseline"="firebrick3", "Cycle1"="gold", "EOT"="dodgerblue4"),
                  tf=c("High"="darkslategray3", "Low"="darkslategray1"),
                  pl=gradient_fun,
                  cb=c("Yes"="steelblue", "No"="tomato"))
  
  onco_t = oncoPrint(data$mat, alter_fun = onco_setting$alter_fun, col = onco_setting$col, column_title=title,
                     heatmap_legend_param=legend_param,
                     remove_empty_rows = remove_empty_rows,
                     pct_side = "right", row_names_side = "left",
                     right_annotation = rowAnnotation(
                     row_barplot = anno_oncoprint_barplot(
                         type=c("amplification", "gain", "loss", "deletion", "nLOH"))),
                     top_annotation = HeatmapAnnotation(
                     tm=data$meta[colnames(data$mat),"timepoint"],
                     tf=data$meta[colnames(data$mat),"tumor_fraction"],
                     pl=data$meta[colnames(data$mat),"absolute_ploidy"],
                     cb=data$meta[colnames(data$mat),"clin_ben"],
                     col=annot_col))

  if(is.na(column_order[1])){
    clin_ben_strs = data$meta$clin_ben
    column_order = c(onco_t@column_order[which(clin_ben_strs[onco_t@column_order] == "No")],
                     onco_t@column_order[which(clin_ben_strs[onco_t@column_order] == "Yes")])
  }
  if(is.na(row_order[1])){
    return( oncoPrint(data$mat, alter_fun = onco_setting$alter_fun, col = onco_setting$col, column_title=title,
                      heatmap_legend_param=legend_param,
                      column_order=column_order,
                      show_column_names=TRUE,
                      remove_empty_rows = remove_empty_rows,
                      pct_side = "right", row_names_side = "left",
                      pct_gp=gpar(fontsize=6),
                      column_names_gp=gpar(fontsize=5),
                      right_annotation = rowAnnotation(
                      row_barplot = anno_oncoprint_barplot(
                         type=c("amplification", "gain", "loss", "deletion", "nLOH"))),
                      top_annotation = HeatmapAnnotation(
                      tm=data$meta[colnames(data$mat),"timepoint"],
                      tf=data$meta[colnames(data$mat),"tumor_fraction"],
                      pl=data$meta[colnames(data$mat),"absolute_ploidy"],
                      cb=data$meta[colnames(data$mat),"clin_ben"],
                      col=annot_col, simple_anno_size = unit(2, "mm"))) )  
  } else {
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
                         type=c("amplification", "gain", "loss", "deletion", "nLOH"))),
                      top_annotation = HeatmapAnnotation(
                      tm=data$meta[colnames(data$mat),"timepoint"],
                      tf=data$meta[colnames(data$mat),"tumor_fraction"],
                      pl=data$meta[colnames(data$mat),"absolute_ploidy"],
                      cb=data$meta[colnames(data$mat),"clin_ben"],
                      col=annot_col, simple_anno_size = unit(2, "mm"))) )
  }
}
