library(data.table)
library(tidyverse)
library(DESeq2)

setwd('C:/Users/Timur/projects/nug')

## data

df <- fread('data/raw_counts.csv')
df <- column_to_rownames(df, 'V1')

colnames(df) <- gsub('DHA200', 'DHA', colnames(df))

df <- df[, c('LUC_Ctrl_1', 'LUC_Ctrl_2', 'LUC_Ctrl_3', 
             'IGFBP6_Ctrl_1', 'IGFBP6_Ctrl_2', 'IGFBP6_Ctrl_3',
             'LUC_DHA_1', 'LUC_DHA_2', 'LUC_DHA_3', 
             'IGFBP6_DHA_1', 'IGFBP6_DHA_2', 'IGFBP6_DHA_3',
             'LUC_Erastin_1', 'LUC_Erastin_2', 'LUC_Erastin_3', 
             'IGFBP6_Erastin_1', 'IGFBP6_Erastin_2')]

samples <- data.frame(Sample = colnames(df),
                      Cell.Line = c(rep('LUC', 3),
                                    rep('IGFBP6', 3),
                                    rep('LUC', 3),
                                    rep('IGFBP6', 3),
                                    rep('LUC', 3),
                                    rep('IGFBP6', 2)),
                      Condition = c(rep('Ctrl', 3),
                                      rep('Ctrl', 3),
                                      rep('DHA', 3),
                                      rep('DHA', 3),
                                      rep('Erastin', 3),
                                      rep('Erastin', 2)))


samples$Cell.Line <- factor(samples$Cell.Line)
samples$Condition <- factor(samples$Condition)

## DEA

run_DEA <- function(cell_line, condition1, condition2, df, samples) {
    subset_samples <- samples %>% filter(Cell.Line == cell_line)
    subset_df <- df[, subset_samples$Sample]
    
    dds <- DESeqDataSetFromMatrix(subset_df, 
                                  colData = subset_samples, 
                                  design = ~ Condition)
    dds <- DESeq(dds)
    
    res <- results(dds, alpha = 0.05, 
                   contrast = c('Condition', condition1, condition2))
    
    resLFC <- lfcShrink(dds, 
                        contrast = c('Condition', condition1, condition2), 
                        res = res, 
                        type = 'normal')
    
    resLFC <- as.data.frame(resLFC)
    resLFC <- rownames_to_column(resLFC, 'gene')
    resLFC <- resLFC[resLFC$padj < 0.05,]
    resLFC <- resLFC[complete.cases(resLFC$log2FoldChange),]
    
    return(resLFC)
}


# res_LUC_DHA200_vs_Ctrl <- run_DEA('LUC', 'DHA200', 'Ctrl', df, samples)
# 
# res_LUC_Erastin_vs_Ctrl <- run_DEA('LUC', 'Erastin', 'Ctrl', df, samples)

res_IGFBP6_DHA200_vs_Ctrl <- run_DEA('IGFBP6', 'DHA', 'Ctrl', df, samples)

res_IGFBP6_Erastin_vs_Ctrl <- run_DEA('IGFBP6', 'Erastin', 'Ctrl', df, samples)

fwrite(as.data.frame(res_LUC_DHA200_vs_Ctrl), 'dea/luc_dha200_vs_ctrl.csv')
fwrite(as.data.frame(res_LUC_Erastin_vs_Ctrl), 'dea/luc_erastin_vs_ctrl.csv')
fwrite(as.data.frame(res_IGFBP6_DHA200_vs_Ctrl), 'dea/IGFBP6_dha200_vs_ctrl.csv')
fwrite(as.data.frame(res_IGFBP6_Erastin_vs_Ctrl), 'dea/IGFBP6_erastin_vs_ctrl.csv')

## GSEA

library(fgsea)
library(ggplot2)
library(org.Hs.eg.db)

# gsea <- function(res, gmt_path) {
#     ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
#                                         key=res$gene, 
#                                         columns='SYMBOL',
#                                         keytype='ENSEMBL')
#     
#     ens2symbol <- as_tibble(ens2symbol)
#     
#     res <- inner_join(res, ens2symbol, by=c('gene'='ENSEMBL'))
#     
#     res2 <- res %>% 
#         dplyr::select(SYMBOL, stat) %>% 
#         na.omit() %>% 
#         distinct() %>% 
#         group_by(SYMBOL) %>% 
#         summarize(stat=mean(stat))
#     
#     ranks <- deframe(res2)
#     
#     pathways <- gmtPathways(gmt_path)
#     
#     fgseaRes <- fgseaMultilevel(pathways = pathways, 
#                                 stats = ranks,
#                                 nproc = 1)
#     
#     fwrite(fgseaRes, paste0('gsea/', df_path))
#     
#     fgseaResTidy <- fgseaRes %>%
#         as_tibble() %>%
#         arrange(desc(NES))
#     
#     jpeg(paste0('gsea/plots/', str_split(df_path, '\\.')[[1]][1], '.jpeg'))
#     ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#         geom_col(aes(fill=padj<0.05)) +
#         coord_flip() +
#         labs(x='Pathway', y='Normalized Enrichment Score',
#              title='Hallmark pathways NES') + 
#         theme_minimal()
#     dev.off()
# }

# dfs <- list.files('dea/')

# for (df_path in dfs) {
#     res <- fread(paste0('dea/', df_path))
#     res$gene <- gsub('\\.\\d+$', '', res$gene)
#     
#     gsea(res, 'data/gmt/h.all.v2024.1.Hs.symbols.gmt')
# }

dfs <- list.files('dea/')

### IGFBP6_dha200_vs_ctrl

res <- fread('dea/IGFBP6_dha200_vs_ctrl.csv')
res$gene <- gsub('\\.\\d+$', '', res$gene)

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$gene, 
                                    columns='SYMBOL',
                                    keytype='ENSEMBL')

ens2symbol <- as_tibble(ens2symbol)

res <- inner_join(res, ens2symbol, by=c('gene'='ENSEMBL'))

res2 <- res %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat))

ranks <- deframe(res2)

pathways <- gmtPathways('data/gmt/h.all.v2024.1.Hs.symbols.gmt')

fgseaRes <- fgseaMultilevel(pathways = pathways, 
                            stats = ranks,
                            nproc = 1)

fwrite(fgseaRes, paste0('gsea/', 'IGFBP6_dha200_vs_ctrl.csv'))

fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x='Pathway', y='Normalized Enrichment Score',
         title='IGFBP6 DHA200 vs Ctrl') + 
    theme_minimal()
ggsave('gsea/plots/IGFBP6_dha200_vs_ctrl.jpeg',
       width = 8, height = 8)

### IGFBP6_erastin_vs_ctrl

res <- fread(paste0('dea/', dfs[2]))
res$gene <- gsub('\\.\\d+$', '', res$gene)

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$gene, 
                                    columns='SYMBOL',
                                    keytype='ENSEMBL')

ens2symbol <- as_tibble(ens2symbol)

res <- inner_join(res, ens2symbol, by=c('gene'='ENSEMBL'))

res2 <- res %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat))

ranks <- deframe(res2)

pathways <- gmtPathways('data/gmt/h.all.v2024.1.Hs.symbols.gmt')

fgseaRes <- fgseaMultilevel(pathways = pathways, 
                            stats = ranks,
                            nproc = 1)

fwrite(fgseaRes, paste0('gsea/', dfs[2]))

fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x='Pathway', y='Normalized Enrichment Score',
         title='IGFBP6 Erastin vs Ctrl') + 
    theme_minimal()
ggsave('gsea/plots/IGFBP6_erastin_vs_ctrl.jpeg',
       width = 8, height = 8)


### luc_dha200_vs_ctrl

res <- fread(paste0('dea/', dfs[3]))
res$gene <- gsub('\\.\\d+$', '', res$gene)

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$gene, 
                                    columns='SYMBOL',
                                    keytype='ENSEMBL')

ens2symbol <- as_tibble(ens2symbol)

res <- inner_join(res, ens2symbol, by=c('gene'='ENSEMBL'))

res2 <- res %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat))

ranks <- deframe(res2)

pathways <- gmtPathways('data/gmt/h.all.v2024.1.Hs.symbols.gmt')

fgseaRes <- fgseaMultilevel(pathways = pathways, 
                            stats = ranks,
                            nproc = 1)

fwrite(fgseaRes, paste0('gsea/', dfs[3]))

fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x='Pathway', y='Normalized Enrichment Score',
         title='LUC DHA200 vs Ctrl') + 
    theme_minimal()
ggsave('gsea/plots/luc_dha200_vs_ctrl.jpeg',
       width = 8, height = 8)

### luc_erastin_vs_ctrl

res <- fread(paste0('dea/', dfs[4]))
res$gene <- gsub('\\.\\d+$', '', res$gene)

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$gene, 
                                    columns='SYMBOL',
                                    keytype='ENSEMBL')

ens2symbol <- as_tibble(ens2symbol)

res <- inner_join(res, ens2symbol, by=c('gene'='ENSEMBL'))

res2 <- res %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat))

ranks <- deframe(res2)

pathways <- gmtPathways('data/gmt/h.all.v2024.1.Hs.symbols.gmt')

fgseaRes <- fgseaMultilevel(pathways = pathways, 
                            stats = ranks,
                            nproc = 1)

fwrite(fgseaRes, paste0('gsea/', dfs[4]))

fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x='Pathway', y='Normalized Enrichment Score',
         title='LUC Erastin vs Ctrl') + 
    theme_minimal()
ggsave('gsea/plots/luc_erastin_vs_ctrl.jpeg',
       width = 8, height = 8)


## volcano plots

library(EnhancedVolcano)

volcano <- function(cell_line, condition1, condition2, df, samples) {
    subset_samples <- samples %>% filter(Cell.Line == cell_line)
    subset_df <- df[, subset_samples$Sample]
    
    dds <- DESeqDataSetFromMatrix(subset_df, 
                                  colData = subset_samples, 
                                  design = ~ Condition)
    dds <- DESeq(dds)
    
    res <- results(dds, contrast = c('Condition', condition1, condition2))
    
    res <- lfcShrink(dds, contrast = c('Condition', condition1, condition2), 
                     res = res, type = 'normal')
    
    res <- as.data.frame(res)

    res$gene <- gsub('\\.\\d+$', '', rownames(res))

    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        key = res$gene,
                                        columns = 'SYMBOL',
                                        keytype = 'ENSEMBL')

    ens2symbol <- as_tibble(ens2symbol)

    res <- inner_join(res, ens2symbol, by = c('gene' = 'ENSEMBL'))
    
    res <- res[!(is.na(res$SYMBOL)),]
    
    res_ <- copy(res)
    
    top_genes_padj <- res_$SYMBOL[order(res_$padj)][1:5]
    
    res_ <- res_[!(res_$SYMBOL %in% top_genes_padj),]
    
    top_genes_lfc_pos <- res_$SYMBOL[order(res_$log2FoldChange, decreasing = TRUE)][1:4]
    
    res_ <- res_[!(res_$SYMBOL %in% top_genes_lfc_pos),]
    
    top_genes_lfc_neg <- res_$SYMBOL[order(res_$log2FoldChange)][1:4]
    
    top_genes <- unique(c(top_genes_padj, top_genes_lfc_pos, top_genes_lfc_neg))
    
    print(top_genes)
    
    p <- EnhancedVolcano(res,
                         lab = res$SYMBOL,
                         selectLab = top_genes,
                         x = 'log2FoldChange',
                         y = 'padj',
                         caption = NULL,
                         drawConnectors = TRUE,
                         typeConnectors = 'open',
                         boxedLabels = F,
                         legendLabels = c("NS", expression(Log[2] ~ FC), "padj", expression(padj ~ and
                                                                                               ~ log[2] ~ FC)),
                         title = paste(paste0('sh-', cell_line), condition1, 'vs', condition2),
                         subtitle = NULL,
                         titleLabSize = 28,
                         legendLabSize = 21,
                         axisLabSize = 24,
                         labSize = 7,
                         FCcutoff = log2(1.5),
                         pCutoff = 0.05,
                         labFace = 'bold')
    
    
    ggsave(paste0('volcano_plots/', 
                  paste(cell_line, condition1, 'vs', condition2, sep = '_'), 
                  '.jpeg'), 
           plot = p, 
           width = 8, height = 8)
}

volcano('LUC', 'DHA', 'Ctrl', df, samples)

volcano('LUC', 'Erastin', 'Ctrl', df, samples)

volcano('IGFBP6', 'DHA', 'Ctrl', df, samples)

volcano('IGFBP6', 'Erastin', 'Ctrl', df, samples)

















library(EnhancedVolcano)

volcano <- function(cell_line, condition1, condition2, df, samples,
                    highlight_genes = NULL) {
    subset_samples <- samples %>% filter(Cell.Line == cell_line)
    subset_df <- df[, subset_samples$Sample]
    
    dds <- DESeqDataSetFromMatrix(subset_df,
                                  colData = subset_samples,
                                  design = ~ Condition)
    dds <- DESeq(dds)
    
    res <- results(dds, contrast = c('Condition', condition1, condition2))
    res <- lfcShrink(dds, contrast = c('Condition', condition1, condition2),
                     res = res, type = 'normal') %>% as.data.frame()
    res$gene <- gsub('\\.\\d+$', '', rownames(res))
    
    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        key = res$gene,
                                        columns = 'SYMBOL',
                                        keytype = 'ENSEMBL') %>% as_tibble()
    
    res <- inner_join(res, ens2symbol, by = c('gene' = 'ENSEMBL')) %>%
        filter(!is.na(SYMBOL))
    
    if (is.null(highlight_genes)) highlight_genes <- character(0)
    
    down_res <- res %>% filter(log2FoldChange < 0)
    
    fc_rank   <- down_res$SYMBOL[order(down_res$log2FoldChange)]
    padj_rank <- down_res$SYMBOL[order(down_res$padj)]
    
    down_highlight <- fc_rank[1:3]                        # топ-3 по FC
    down_highlight <- c(down_highlight,                  # + топ-2 по padj,
                        setdiff(padj_rank, down_highlight)[1:2])
    
    if (!'SCD' %in% down_highlight)                       # + SCD
        down_highlight <- c(down_highlight, 'SCD')
    
    # если дубли всё-таки сократили список < 6, добираем следующими генами
    if (length(down_highlight) < 6) {
        extras <- setdiff(c(fc_rank, padj_rank), down_highlight)
        down_highlight <- c(down_highlight, extras[1:(6 - length(down_highlight))])
    }
    
    
    up_highlight <- intersect(highlight_genes,
                              res$SYMBOL[res$log2FoldChange > 0])
    
    select_genes <- unique(c(down_highlight, up_highlight))
    select_genes <- select_genes[select_genes %in% res$SYMBOL]
    
    p <- EnhancedVolcano(res,
                         lab = res$SYMBOL,
                         selectLab = select_genes,
                         x = 'log2FoldChange',
                         y = 'padj',
                         caption = NULL,
                         drawConnectors = TRUE,
                         typeConnectors = 'open',
                         boxedLabels = TRUE,
                         legendLabels = c('NS', expression(Log[2] ~ FC), 'padj',
                                          expression(padj ~ and ~ log[2] ~ FC)),
                         title = paste(paste0('sh-', cell_line), condition1, 'vs', condition2),
                         subtitle = NULL,
                         titleLabSize = 28,
                         legendLabSize = 21,
                         axisLabSize = 24,
                         # FCcutoff = log2(1.5),
                         labSize = 7,
                         FCcutoff = log2(1.5),
                         pCutoff = 0.05,
                         labFace = 'bold',
                         cutoffLineType = 'twodash',
                         legendPosition = 'right',
                         col = c('grey65', 'blue', 'darkorange1', 'forestgreen'),
                         colAlpha = 0.9,
                         gridlines.minor = FALSE)
    
    ggsave(paste0('volcano_plots/new/',
                  paste(cell_line, condition1, 'vs', condition2, sep = '_'),
                  '.jpeg'),
           plot = p,
           width = 12, height = 8)
}



genes_to_bold <- c('HMOX1', 'SLC7A11', 'NQO1', 'FTH1', 'GCLM', 'FTL',
                   'SPTSSB', 'LIPG', 'DMBT1', 'OR2M3', 'PDE4B', 'SCD')

volcano('IGFBP6', 'DHA', 'Ctrl', df, samples,
        highlight_genes = genes_to_bold)

volcano('IGFBP6', 'Erastin', 'Ctrl', df, samples,
        highlight_genes = genes_to_bold)



library(EnhancedVolcano)

volcano <- function(cell_line, condition1, condition2, df, samples,
                    highlight_genes = NULL,
                    n_fc = 3, n_padj = 3,
                    fc_cutoff = log2(1.5), p_cutoff = 0.05,
                    save_dir = 'volcano_plots/new', width = 12, height = 8) {
    
    subset_samples <- samples %>% filter(Cell.Line == cell_line)
    subset_df      <- df[, subset_samples$Sample]
    
    dds <- DESeqDataSetFromMatrix(subset_df,
                                  colData = subset_samples,
                                  design = ~ Condition) %>% DESeq()
    
    res <- results(dds, contrast = c('Condition', condition1, condition2))
    res <- lfcShrink(dds, contrast = c('Condition', condition1, condition2),
                     res = res, type = 'normal') %>% as.data.frame()
    res$gene <- gsub('\\.\\d+$', '', rownames(res))
    
    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        key     = res$gene,
                                        columns = 'SYMBOL',
                                        keytype = 'ENSEMBL') %>% as_tibble()
    
    res <- res %>% inner_join(ens2symbol, by = c('gene' = 'ENSEMBL')) %>%
        filter(!is.na(SYMBOL))
    
    get_top <- function(tbl, desc_fc) {
        fc_order   <- if (desc_fc) order(-tbl$log2FoldChange) else order(tbl$log2FoldChange)
        fc_genes   <- tbl$SYMBOL[fc_order][1:n_fc]
        padj_genes <- tbl$SYMBOL[order(tbl$padj)][1:n_padj]
        top <- unique(c(fc_genes, padj_genes))
        if (length(top) < n_fc + n_padj) {
            extras <- setdiff(tbl$SYMBOL[order(tbl$padj)], top)
            top <- c(top, extras[1:((n_fc + n_padj) - length(top))])
        }
        top
    }
    
    up_high   <- get_top(res %>% filter(log2FoldChange > 0),  TRUE)
    down_high <- get_top(res %>% filter(log2FoldChange < 0), FALSE)
    
    select_genes <- unique(c(up_high, down_high,
                             if (!is.null(highlight_genes)) highlight_genes else character(0)))
    select_genes <- select_genes[select_genes %in% res$SYMBOL]
    
    p <- EnhancedVolcano(res,
                         lab            = res$SYMBOL,
                         selectLab      = select_genes,
                         x              = 'log2FoldChange',
                         y              = 'padj',
                         caption        = NULL,
                         drawConnectors = TRUE,
                         typeConnectors = 'open',
                         boxedLabels    = TRUE,
                         legendLabels   = c('NS', expression(Log[2]~FC), 'padj',
                                            expression(padj~and~log[2]~FC)),
                         title          = paste(paste0('sh-', cell_line),
                                                condition1, 'vs', condition2),
                         titleLabSize   = 28,
                         legendLabSize  = 21,
                         subtitle = NULL,
                         axisLabSize    = 24,
                         labSize        = 7,
                         FCcutoff       = fc_cutoff,
                         pCutoff        = p_cutoff,
                         labFace        = 'bold',
                         cutoffLineType = 'twodash',
                         legendPosition = 'right',
                         col            = c('grey65', 'blue', 'darkorange1', 'forestgreen'),
                         colAlpha       = 0.9,
                         gridlines.minor = FALSE)
    
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(file.path(save_dir,
                     paste0(paste(cell_line, condition1, 'vs', condition2, sep = '_'),
                            '.jpeg')),
           plot  = p,
           width = width, height = height)
}





volcano('IGFBP6', 'DHA', 'Ctrl', df, samples)

volcano('IGFBP6', 'Erastin', 'Ctrl', df, samples)

























library(data.table)
library(dplyr)
library(ggvenn)
library(ggplot2)

dir.create('venn', showWarnings = FALSE)

conds  <- c('dha200', 'erastin')
pretty <- c('ДГК 200 мкМ', 'Эрастин')
names(pretty) <- conds

files <- list(
    dha200  = 'dea/IGFBP6_dha200_prep.csv',
    erastin = 'dea/IGFBP6_erastin_prep.csv'
)

up_list   <- list()
down_list <- list()

for (cnd in conds) {
    tbl   <- fread(files[[cnd]], data.table = FALSE)
    label <- pretty[[cnd]]
    up_list  [[label]] <- tbl %>% filter(FC > 0) %>% pull(symbol) %>% unique()
    down_list[[label]] <- tbl %>% filter(FC < 0) %>% pull(symbol) %>% unique()
}

make_logical_df <- function(lst) {
    allg <- unique(unlist(lst))
    data.frame(lapply(lst, \(v) allg %in% v), row.names = allg, check.names = FALSE)
}

cols <- c('#4E79A7', '#E15759')

build_venn <- function(lst, title, file) {
    p <- ggvenn(
        data          = make_logical_df(lst),
        fill_color    = cols,
        stroke_size   = .5,
        set_name_size = 8,
        text_size     = 8
    ) +
        ggtitle(title) +
        theme(
            plot.title = element_text(size = 26, face = 'bold',
                                      hjust = .5, vjust = -4)
        )
    ggsave(file, p, width = 10, height = 8, dpi = 300)
}

build_venn(up_list,   'Гены с повышенной экспрессией',   'venn/up_genes.jpeg')
build_venn(down_list, 'Гены с пониженной экспрессией', 'venn/down_genes.jpeg')











library(dplyr)
library(stringr)
library(ggvenn)
library(ggplot2)

dir.create('venn_pathways', showWarnings = FALSE)

conds  <- c('dha200', 'erastin')
pretty <- c('ДГК 200 мкМ', 'Эрастин'); names(pretty) <- conds
files  <- list(
    dha200  = 'gsea/igfbp6_dha200_vs_ctrl.csv',
    erastin = 'gsea/igfbp6_erastin_vs_ctrl.csv'
)

q_cut <- .25            # порог FDR

read_clean <- \(f) {
    fread(f, data.table = FALSE) |>
        mutate(
            across(c(padj, NES),
                   \(x) as.numeric(str_replace_all(x, ',', '.'))  # ',' → '.'
            ),
            pathway = str_trim(pathway)                      # убираем пробелы
        )
}

up   <- list()
down <- list()

for (cnd in conds) {
    tbl <- read_clean(files[[cnd]])
    lab <- pretty[[cnd]]
    up  [[lab]] <- tbl |> filter(padj < q_cut,  NES > 0) |> pull(pathway) |> unique()
    down[[lab]] <- tbl |> filter(padj < q_cut,  NES < 0) |> pull(pathway) |> unique()
}

message('↑ DHA: ', length(up[[pretty['dha200']] ]),
        ', ↑ Erastin: ', length(up[[pretty['erastin']]]))

make_log <- \(lst) {
    all <- unique(unlist(lst))
    data.frame(lapply(lst, \(v) all %in% v), row.names = all,
               check.names = FALSE)
}

cols <- c('#4E79A7', '#E15759')

save_venn <- \(lst, title, file) {
    p <- ggvenn(
        data          = make_log(lst),
        fill_color    = cols,
        stroke_size   = .5,
        set_name_size = 8,
        text_size     = 8
    ) +
        ggtitle(title) +
        theme(plot.title = element_text(size = 26, face = 'bold',
                                        hjust = .5, vjust = -4))
    
    ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
}


save_venn(up,   'Активированные пути', 'venn_pathways/up_pathways.jpeg')
save_venn(down, 'Подавленные пути', 'venn_pathways/down_pathways.jpeg')

