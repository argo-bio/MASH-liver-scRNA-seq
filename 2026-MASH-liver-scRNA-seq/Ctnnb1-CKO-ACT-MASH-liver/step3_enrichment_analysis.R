## May 6 2026

### step3 CKO vs WT and ACT vs WT DEG analysis and enrichment analysis

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(readxl)
deg_CKO_vs_WT <- read.csv('./output/deg/group_CKO-vs-WT-all_diffexp_genes_anno.csv')
deg_ACT_vs_WT <- read.csv('./output/deg/group_ACT-vs-WT-all_diffexp_genes_anno.csv')

logfc_cutoff <- 0.25
padj_cutoff <- 0.05

CKO_up <- deg_CKO_vs_WT %>%
  filter(log2FoldChange > logfc_cutoff, p.value < padj_cutoff, gene_type=='protein_coding') %>%
  pull(gene) %>%
  unique()

CKO_down <- deg_CKO_vs_WT %>%
  filter(log2FoldChange < -logfc_cutoff, p.value < padj_cutoff, gene_type=='protein_coding') %>%
  pull(gene) %>%
  unique()

ACT_up <- deg_ACT_vs_WT %>%
  filter(log2FoldChange > logfc_cutoff, p.value < padj_cutoff, gene_type=='protein_coding') %>%
  pull(gene) %>%
  unique()

ACT_down <- deg_ACT_vs_WT %>%
  filter(log2FoldChange < -logfc_cutoff, p.value < padj_cutoff, gene_type=='protein_coding') %>%
  pull(gene) %>%
  unique()

write.csv(
  data.frame(gene = CKO_up_ACT_down),
  "./output/deg/CKO_up_ACT_down_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = ACT_up_CKO_down),
  "./output/deg/ACT_up_CKO_down_genes.csv",
  row.names = FALSE
)

ACT_up_CKO_down <- read.csv('./output/deg/ACT_up_CKO_down_genes_filter.csv')
ACT_up_CKO_down <- ACT_up_CKO_down$gene

### enrichment analysis
gene_df <- bitr(
  ACT_up_CKO_down,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

entrez_genes <- unique(gene_df$ENTREZID)

#### GO enrichment
ego_bp <- enrichGO(
  gene          = entrez_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

write.csv(
  as.data.frame(ego_bp),
  "ACT_up_CKO_down_GO_BP_enrichment.csv",
  row.names = FALSE
)

#### KEGG enrichment
ekegg <- enrichKEGG(
  gene          = entrez_genes,
  organism      = "mmu",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

ekegg <- setReadable(
  ekegg,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID"
)

write.csv(
  as.data.frame(ekegg),
  "ACT_up_CKO_down_KEGG_enrichment.csv",
  row.names = FALSE
)

pdf("ACT_up_CKO_down_GO_BP_dotplot.pdf", width = 8, height = 6)
dotplot(ego_bp, showCategory = 20) +
  ggtitle("GO BP enrichment: ACT up & CKO down genes")
dev.off()

ego_df <- as.data.frame(ego_bp)
write.csv(ego_df,file = './output/deg/act_up_cko_down_go_enrichment.csv')

ego_df2 <- read.csv('./output/deg/act_up_cko_down_go_enrichment_filter.csv',row.names = 'X')
ego_df2$logp <- -log10(ego_df2$p.adjust)


ggplot(ego_df2,aes(y=reorder(Description, zScore),x=RichFactor))+
  
  geom_point(aes(size=Count,color=logp))+
  
  scale_color_gradient(low = "red",high ="blue")+
  scale_size_continuous(range = c(5, 8))+
  
  labs(color=expression(logp,size="Count"),
       
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()
ggsave('./output/deg/act_up_cko_down_go_enrichment_filter.pdf',width = 8,height = 6)

pdf("ACT_up_CKO_down_KEGG_dotplot.pdf", width = 8, height = 6)
dotplot(ekegg, showCategory = 20) +
  ggtitle("KEGG enrichment: ACT up & CKO down genes")
dev.off()

summary_table <- bind_rows(
  data.frame(Gene = ACT_up,    Group = "ACT vs WT", Regulation = "Upregulated"),
  data.frame(Gene = ACT_down,  Group = "ACT vs WT", Regulation = "Downregulated"),
  data.frame(Gene = CKO_up,    Group = "CKO vs WT", Regulation = "Upregulated"),
  data.frame(Gene = CKO_down,  Group = "CKO vs WT", Regulation = "Downregulated"))

write.csv(summary_table,'./output/deg/DEGs_ACTvsWT_CKOvsWT.csv')


### vocanol plot
library(dplyr)
library(ggplot2)
library(ggrepel)

rm(list = ls())

deg_CKO_vs_WT <- read.csv('./output/deg/group_CKO-vs-WT-all_diffexp_genes.csv')
deg_ACT_vs_WT <- read.csv('./output/deg/group_ACT-vs-WT-all_diffexp_genes.csv')

logfc_cutoff <- 0.25
padj_cutoff <- 0.05

plot_df <- deg_CKO_vs_WT %>%
  dplyr::select(
    gene,
    CKO_log2FC = log2FoldChange,
    CKO_padj = p.value
  ) %>%
  inner_join(
    deg_ACT_vs_WT %>%
      dplyr::select(
        gene,
        ACT_log2FC = log2FoldChange,
        ACT_padj = p.value
      ),
    by = "gene"
  )


label_genes <- c(
  "Scd1", "Trem2","Srebf1","Csnk1a1","Wnk1","Fasn","Acaca"
)

plot_df$label <- ifelse(
  plot_df$gene %in% label_genes &
    plot_df$group != "Other",
  plot_df$gene,
  NA
)

plot_df <- plot_df %>% filter(abs(CKO_log2FC) < 2, CKO_padj<0.05,
                              abs(ACT_log2FC) < 2, ACT_padj<0.05)

ggplot(plot_df,
       aes(x = ACT_log2FC,
           y = CKO_log2FC)) +
  geom_point(aes(color = group),size = 1.8,alpha = 0.7) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_hline(
    yintercept = c(-logfc_cutoff, logfc_cutoff),
    linetype = "dotted",
    color = "grey70"
  ) +
  geom_vline(
    xintercept = c(-logfc_cutoff, logfc_cutoff),
    linetype = "dotted",
    color = "grey70"
  ) +
  geom_text_repel(
    aes(label = label),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color='red',
    nudge_x = 0.5,
    nudge_y = -0.5,
    force=5,
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  scale_color_manual(
    values = c(
      "CKO up & ACT down" = "#D73027",
      "ACT up & CKO down" = "#4575B4",
      "Other" = "grey80"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "",
    x = "log2FC (ACT vs WT)",
    y = "log2FC (CKO vs WT)",
    color = NULL
  ) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
    )
ggsave('./output/deg/dual_vovannol.pdf',width = 8,height = 5)
