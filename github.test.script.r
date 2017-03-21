###############################################################################
# Make heatmap function                                                       #
###############################################################################

make.hm <- function(
  m.df1, 
  filename = "heatmap", 
  k.number = 10, 
  n.colors = 1000, 
  hclust.method = "complete", 
  dist.method = "euclidean", 
  main = "") {
  library(RColorBrewer)
  library(gplots)
  hclustfunc <- function(x) hclust(x, method = hclust.method)
  distfunc <- function(x) dist(x, method = dist.method)
  d <- distfunc(m.df1)
  fit <- hclustfunc(d)
  clusters <- cutree(fit, k = k.number)
  nofclust.height <- length(unique(as.vector(clusters)))
  a = max(abs(m.df1))
  breaks = seq((-1 * a), a, length.out = n.colors)
  hmcols <- colorRampPalette(c("blue", "white", "red"))(n = (length(breaks) - 
                                                               1))
  selcol <- colorRampPalette(brewer.pal(12, "Set3"))
  selcol2 <- colorRampPalette(brewer.pal(9, "Set1"))
  clustcol.height = selcol2(nofclust.height)
  pdf(paste(filename, "pdf", sep = "."))
  hm = heatmap.2(
    m.df1, 
    trace = "none", 
    dendrogram = "both", 
    density.info = "none", 
    keysize = 1, 
    key = TRUE, 
    Colv = F, 
    hclust = hclustfunc, 
    distfun = distfunc, 
    col = hmcols, 
    symbreak = T, 
    labRow = rep("", nrow(m.df1)),
    labCol = gsub("_", " ",gsub("lg2_avg_", "", colnames(m.df1))),
    RowSideColors = clustcol.height[clusters], 
    margins = c(10, 10), cexCol = 1, cexRow = 0.5, srtCol = 45, 
    srtRow = 0, 
    main = main, 
    breaks = breaks, 
    sepcolor = "black", 
    sepwidth = c(5e-04, 5e-05), 
    colsep = c(0, ncol(m.df1)), 
    rowsep = c(0, nrow(m.df1)))
  dev.off()
  sorted = m.df1[match(rev(labels(hm$rowDendrogram)), rownames(m.df1)), 
                 ]
  sorted = sorted[, hm$colInd]
  pdf(paste(filename, "colorkey.pdf", sep = "."))
  plot.new()
  par(lend = 1)
  legend("topleft", legend = 1:nofclust.height, col = clustcol.height, 
         lty = 1, lwd = 10)
  dev.off()
  df.res = list(sorted = sorted, clusters = clusters)
  return(df.res)
  return(clusters)
}

df.cluster[df.cluster > hm.cut.off] = hm.cut.off
df.cluster[df.cluster < (-1) * hm.cut.off] = (-1) * hm.cut.off
m.cluster = data.matrix(df.cluster)
m.cluster[is.na(m.cluster)] = 0
#m.cluster[m.cluster == 0] = 0.01
hm.res = make.hm(m.cluster, 
                 filename = "heatmap", 
                 k.number = n.hm.cluster, 
                 n.colors = 1000, 
                 hclust.method = "ward.D2", 
                 dist.method = "euclidean", 
                 main = "")

#if (!use.logFC.columns.for.heatmap) {
#  new.sample.order.lg2.avg = colnames(hm.res$sorted)
#  new.sample.order.norm.counts = gsub("lg2_avg_", "norm_counts_", 
#                                      colnames(hm.res$sorted))
#  other.cols = names(df.data)[-match(new.sample.order.norm.counts, 
#                                     names(df.data))]
#  other.cols = other.cols[-match(new.sample.order.lg2.avg, 
#                                 other.cols)]
#  new.order = append(other.cols, new.sample.order.norm.counts)
#  new.order = append(new.order, new.sample.order.lg2.avg)
#  if (sum(new.order %in% names(df.data)) == length(names(df.data))) {
#    df.data = df.data[, new.order]
#  }
#}

df.clust.order = data.frame(hm.res$sorted)
cluster.ordered.sample.vector = names(df.clust.order)
df.clust.order[["cluster_order"]] = 1:nrow(df.clust.order)
df.clust.order[, "row_id"] = row.names(df.clust.order)
df.clust.order = df.clust.order[, c("row_id", "cluster_order")]
remove = as.vector(na.omit(match(df.clust.order[, "row_id"], 
                                 df.data[, "row_id"])))
id.vector = as.vector(df.data[-remove, "row_id"])
df.rest = data.frame(id.vector, rep(0, length(id.vector)))
names(df.rest) = names(df.clust.order)
df.clust.order = rbind(df.clust.order, df.rest)
df.data = merge(df.data, df.clust.order, by.x = "row_id", 
                by.y = "row_id")
df.data = df.data[!is.na(df.data[, gene.id.column]), ]
df.data = unique(df.data)

# Add cluster id
df.cluster.id <- data.frame(na.omit(hm.res$clusters))
df.cluster.id[["row_id"]] <- row.names(df.cluster.id)
names(df.cluster.id)<- c("cluster_id", "row_id")
df.cluster.id <- df.cluster.id[grep("R", df.cluster.id$row_id),]
# Adding all other ids
row_id <- df.data[!(df.data$row_id %in% df.cluster.id$row_id), "row_id"]
cluster_id <- rep(0, length(row_id))
df.add <- rbind(df.cluster.id, data.frame(cluster_id, row_id))

df.data <- merge(df.data, df.add, by.x = "row_id", by.y="row_id", all=TRUE)
df.data[is.na(df.data)] = ""

# Add gene descripton
df.anno <- read.delim(
  gene.id.table,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Remove all entries from df.anno that are not present in df.data
df.anno <- df.anno[df.anno[,gene.id.column] %in% df.data[,gene.id.column],]
if (!add.uniprot.column){
  df.anno$uniprot = NULL
  df.anno <- unique(df.anno)
}


df.data <- merge(df.data, df.anno, by.x = gene.id.column, by.y = gene.id.column, all=TRUE)

df.data[is.na(df.data)] = ""
df.data = unique(df.data)

if (gene.id.column == "mgi_symbol" | gene.id.column == "ENSMUSG") {
  ENSG <- "ENSMUSG"
} else if (gene.id.column == "hgnc_symbol" | gene.id.column == "ENSG") {
  ENSG <- "ENSG"
}

if (gene.id.column != "mgi_symbol" | gene.id.column != "hgnc_symbol"){
  df.data$gene_description <- paste0(
    df.data$gene_description,
    " (",
    df.data[,gene.id.column],
    ")"
  )
}



if (length(grep("p_site_env", names(df.data))) > 0) {
  #Trim if the sequence window is to big
  length <- nchar(df.data$p_site_env)
  center <- ((length -1)/2)
  df.data$p_site_env <- ifelse(
    (length > 15),
    substr(df.data$p_site_env, center-6,center+8),
    df.data$p_site_env
  )
  
  one = tolower(substr(df.data$p_site_env, 1, 7))
  two = toupper(substr(df.data$p_site_env, 8, 8))
  three = tolower(substr(df.data$p_site_env, 9, 16))
  df.data$p_site_env = paste(one, two, three, sep = "")
  
  ################################################################################
  #Add ppos columns to datatable
  ################################################################################
  ppos.vec = c("ppos_minus_7","ppos_minus_6","ppos_minus_5","ppos_minus_4","ppos_minus_3","ppos_minus_2","ppos_minus_1","ppos",
               "ppos_plus_1", "ppos_plus_2","ppos_plus_3","ppos_plus_4","ppos_plus_5","ppos_plus_6","ppos_plus_7")
  
  
  #In this dataset not all sequences are associated with an p_site_env
  #df.data[df.data$p_site_env == "", "p_site_env"] = substr(df.data[df.data$p_site_env == "", "sequence_window"], 9,23)
  
  for (i in 1:length(ppos.vec))
  {df.data[[ppos.vec[i]]] = sapply(df.data$p_site_env, function(x) substr(x, i,i))
  }
  
  
  # Done adding ppos columns
}
df.data[is.na(df.data)] = ""
df.data = unique(df.data)
df.data = df.data[!is.na(df.data[, gene.id.column]), ]
df.data[["row_names"]] = 1:nrow(df.data)
names(df.data) = gsub("[.]", "_", names(df.data))
names(df.data) = gsub(" ", "_", names(df.data))

return(df.data)
}
## End of function                                                           ##
###############################################################################