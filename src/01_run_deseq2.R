# Importamos paquetes usando library
library (DESeq2)
library (clusterProfiler)
library (ReactomePA)
library (org.Hs.eg.db)
library (enrichplot)
library (ggplot2)


dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

counts_path <- "data/raw/counts.csv.gz" #guardamos string con la ruta relativa del proyecto a donde hemos guardado los datos. Los conteos son los metadatos
metadata_path <- "data/raw/metadata.csv.gz" 
condition_column <- "site"
plate_column <- "plate" # para corregir por placa
condition_test <- "Primary_Solid_Tumor"
condition_ref <- "Solid_Tissue_Normal"
min_total_count <- 20
alpha <- 0.05 # umbral de significancia 
lfc_threshold <- 1 # umbral biológico de log2FoldChange 

message("[1/8] Leyendo datos de entrada (TCGA-BRCA)...")
counts_df <- read.csv(counts_path, row.names = 1, check.names = FALSE) # check.names = FALSE para que no modifique los nombres de las columnas, que son los nombres de las muestras
metadata <- read.csv(metadata_path, row.names = 1, stringsAsFactors = FALSE) # stringsAsFactors = FALSE para que no convierta las columnas de texto en factores. R por defecto convierte los strings en categorías

stopifnot(condition_column %in% colnames(metadata)) # si no se cumple condicion, para aqui. Compara que condition_column site sea una columna
stopifnot(plate_column %in% colnames(metadata)) # plate_column va a marcar los batches 

counts_mat <- as.matrix(counts_df)

message("[2/8] Alineando muestras entre conteos y metadata...")
# quitamos la intersección porque no lo vieron en BMS
#common_samples <- intersect(colnames(counts_mat), rownames(metadata))
#if (length(common_samples) < 2) {
#    stop("No hay suficientes muestras en comun entre counts y metadata.")
#}
metadata <- metadata[colnames(counts_mat), ]

metadata[[condition_column]] <- factor(metadata[[condition_column]])
metadata[[plate_column]] <- factor(metadata[[plate_column]])
metadata[[condition_column]] <- relevel(metadata[[condition_column]], ref = condition_ref)

message("[3/8] Creando objeto DESq2 y filtrando genes de baja cuenta...")
dds <- DESeqDataSetFromMatrix (
    countData = counts_mat,
    colData = metadata, # le pasamos los datos
    design = as.formula(paste("~", condition_column, "+", plate_column)) # ~ significa trata un gen para cada condición (c) y el efecto batch (plate)
)

# Gen con >= min_count lecturas en >= min_samples muestras (tamaño del grupo más pequeño).
# dds <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ] # lo descomentamos para ver las diferencias entre los dos

message("[4/8] Ejecutando DESeq2...")
dds <- DESeq(dds)

# Dispersion: negro = estimacion por gen, rojo = modelo ajustado , azul = valor final.
# Bien si: la nube sigue la linea roja, baja al aumentar los conteos y se estabiliza.
# Mal si: puntos sueltos abajo o la linea se dispara a la izquierda -> filtrar mas. 
png("results/plot_dispersion.png", width = 800, height = 600)
plotDispEsts(dds)
dev.off()

message("[5/8] Extrayendo resultados y guardando tablas...")
res <- results(
    dds,
    contrast = c(condition_column, condition_test, condition_ref),
    alpha = alpha
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.table(
    res_df,
    file = "results/deseq2_all_results.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

norm_df <- data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE)

write.table(
    norm_df,
    file = "results/normalized_counts.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

sig_df <- subset(res_df, !is.na(padj) & padj < alpha)
de_df <- subset(res_df, !is.na(padj) & padj < alpha & abs(log2FoldChange) > lfc_threshold)
up_df <- subset(de_df, log2FoldChange > 0)
down_df <- subset(de_df, log2FoldChange < 0)

top10_up <- head(up_df, 10)
top10_down <- head(down_df, 10)

de_df$top10_label <- "other"
de_df$top10_label[de_df$gene_id %in% top10_up$gene_id] <- "top10_up"
de_df$top10_label[de_df$gene_id %in% top10_down$gene_id] <- "top10_down"

write.table(
    sig_df,
    file = "results/deseq2_sig_padj_0.05.tsv",
    sep = "\t", 
    quote = FALSE,
    row.names = FALSE
)

write.table(
    de_df,
    file = "results/deseq2_de_padj_0.05_lfc_1.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE 
)

write.table(
    top10_up,
    file = "results/top10_up_genes.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    top10_down,
    file = "results/top10_down_genes.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

top10_df <- rbind(
    data.frame(gene_id = top10_up$gene_id, log2FoldChange = top10_up$log2FoldChange,
               direction = "Up (Tumor > Normal)"),
    data.frame(gene_id = top10_down$gene_id, log2FoldChange = top10_down$log2FoldChange,
               direction = "Down (Tumor < Normal)")
)
top10_df$direction <- factor(top10_df$direction,
                             levels = c("Up (Tumor > Normal)", "Down (Tumor < Normal)"))

p_top10 <- ggplot(top10_df, aes(x = reorder(gene_id, log2FoldChange),
                                y = log2FoldChange, fill = direction)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ direction, scales = "free") +
    scale_fill_manual(values = c("Up (Tumor > Normal)" = "#A50F2D",
                                 "Down (Tumor < Normal)" = "#1F5F8B")) +
    labs(x = NULL, y = "log2FoldChange") +
    theme_minimal() +
    theme(legend.position = "none")

ggsave("results/plot_top10_genes.png", p_top10, width = 7, height = 6)

# volcano: genes arriba a los lados = cambio grande y muy significativo.
# Rojo = sube en tumor, azul = baja en tumor, gris = sin cambio. 
volcano_df <- res_df[!is.na(res_df$padj), ]
volcano_df$status <- "No significativo"
volcano_df$status[volcano_df$gene_id %in% up_df$gene_id] <- "Up"
volcano_df$status[volcano_df$gene_id %in% down_df$gene_id] <- "Down"

p_volcano <- ggplot(volcano_df, aes(x= log2FoldChange, y = -log10(padj), color = status)) + 
    geom_point(size = 1) + 
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") + 
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    geom_text(data = head(volcano_df, 10), aes(label = gene_id), vjust = -0.5, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "red", Down = "blue", "No significativo" = "grey")) +
    theme_bw()

ggsave("results/plot_volcano.png", p_volcano, width = 7, height = 6)

# MA: cada punto es un gen; eje x = expresion media, eje y = log2FoldChange.
# La nube debe estar centrada en 0; los puntos azules son significativos. 
png("results/plot_ma.png", width = 800, height = 600)
plotMA(res, ylim = c(-4, 4))
dev.off()
#L: Nube de expresión media centrada en cero. Es un comando de deseq2 donde usamos R básico. 

message("[6/8] Contextualizacion biologica: gen destacado en el contraste...")

top_gene <- de_df[1, ]
writeLines(
    c(
        paste("gene_id:", top_gene$gene_id),
        paste("log2FoldChange:", round(top_gene$log2FoldChange, 4)),
        paste("padj:", signif(top_gene$padj, 3))
    ),
    "results/top_gene_context.txt"
)

message("[7/8] Enriquecimiento GO y Reactome...")
symbol_to_entrez <- bitr(
    de_df$gene_id,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)
universe_entrez <- unique(symbol_to_entrez$ENTREZID)

up_map <- bitr(unique(up_df$gene_id), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_map <- bitr(unique(down_df$gene_id), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Dotplot: procesos/rutas enriquecidos. Mas a la derecha = mas genes del proceso
# color = padj (mas significativo); tamano = numero de genes.

go_up <- enrichGO(
    gene = unique(up_map$ENTREZID),
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
)

write.table(as.data.frame(go_up), "results/go_bp_up.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_go_bp_up.png", dotplot(go_up, showCategory = 10), width = 8, height = 6)
# con showCategory = 10 mostramos los 10 procesos más significativos (solo 10 categorías)

go_down <- enrichGO(
    gene = unique(down_map$ENTREZID),
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
) # hacemos lo mismo para las funciones sobrerepresentadas en los genes que bajan en tumor (infrarrepresentadas)
write.table(as.data.frame(go_down), "results/go_bp_down.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_go_bp_down.png", dotplot(go_down, showCategory = 10), width = 8, height = 6)
# representamos los 10 o 20 primeros y lo vamos ajustando para que no solapen

reactome_up <- enrichPathway(
    gene = unique(up_map$ENTREZID),
    universe = universe_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20,
    readable = TRUE
) # para tener los pathway. 
write.table(as.data.frame(reactome_up), "results/reactome_up.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_reactome_up.png", dotplot(reactome_up, showCategory = 10), width = 8, height = 6)

reactome_down <- enrichPathway(
    gene = unique(down_map$ENTREZID),
    universe = universe_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
)
write.table(as.data.frame(reactome_down), "results/reactome_down.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_reactome_down.png", dotplot(reactome_down, showCategory = 10), width = 8, height = 6)

message("[8/8] Guardando trazabilidad de sesion...")
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")

summary_lines <- c(
    paste("n_samples:", ncol(counts_mat)),
    paste("n_genes_input:", nrow(counts_mat)),
    paste("n_genes_after_filter:", nrow(dds)),
    paste("n_significant_padj_lt", alpha, ":", nrow(sig_df), sep = ""),
    paste("n_de_absLog2FC_gt_", lfc_threshold, "_padj_lt", alpha, ":", nrow(de_df), sep = ""),
    paste("n_up:", nrow(up_df)),
    paste("n_down:", nrow(down_df)),
    paste("top10_up_file:", "results/top10_up_genes.tsv"),
    paste("top10_down_file:", "results/top10_down_genes.tsv"),
    paste("top10_plot_file:", "results/plot_top10_genes.png"),
    paste("volcano_plot_file:", "results/plot_volcano.png")
)

writeLines(summary_lines, "logs/run_summary.txt")

message("Analisis completado")
