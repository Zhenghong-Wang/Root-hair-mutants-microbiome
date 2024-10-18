#usage:
#otu<-read.csv("HM_otutab_20221124_F1128.csv",header = T,row.names = 1)
#meta<-read.csv("HM_metadata_20221124_F1128.csv",header=T,row.names = 1)
#设置unique(meta$group)的分组的level
#title=""  #设置图的标题
#text_name=FALSE
#otu<-otutabi
#meta<-metai

groupPcoa <- function(otu, meta, betaGroup = "betaGroup", cg = "betaGroup", sg = "betaGroup", title = "title", text_name = FALSE, ellipse = TRUE, formula = TRUE, factor = c("group1", "group2"), point_size = 3) {
  
  colnames(meta) <- ifelse(colnames(meta) == betaGroup, "betaGroup", colnames(meta))
  library(vegan)
  
  idx <- rownames(otu)
  species <- otu[idx, ]
  idx <- which(rowSums(species) > 0)
  species <- na.omit(species[idx,])
  
  min <- min(colSums(species))
  set.seed(315)
  otu <- vegan::rrarefy(t(species), min)
  
  set.seed(221)
  bray_curtis <- vegdist(otu, method = "bray")
  bray <- as.matrix(bray_curtis)
  
  dis <- bray
  dis <- as.dist(dis)
  set.seed(110)
  adonis_res <- adonis2(dis ~ betaGroup, meta, permutations = 999)
  
  library(pairwiseAdonis)
  pairwise.adonis <- pairwise.adonis(x = otu, factors = meta$betaGroup, sim.function = "vegdist", sim.method = "bray", reduce = NULL, perm = 999)
  tab2 <- ggtexttable(pairwise.adonis[, c("pairs", "R2", "p.value", "p.adjusted")], rows = NULL, theme = ttheme("blank")) %>% 
    tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
    tab_add_hline(at.row = nrow(pairwise.adonis) + 1, row.side = "bottom", linewidth = 1)
  
  beta <- bray
  pcoa <- wcmdscale(beta, k = 3, eig = TRUE)
  points <- as.data.frame(pcoa$points) 
  eig <- pcoa$eig
  
  points <- points[order(rownames(points)),]
  meta <- meta[order(rownames(meta)),]
  
  points <- cbind(points, meta[, c("betaGroup", cg, sg)])
  colnames(points) <- c("PCoA1", "PCoA2", "PCoA3", "betaGroup", "cg", "sg") 
  
  library(ggplot2)
  mytheme <- theme(
    plot.title = element_text(color = "black", size = 8, hjust = 0.5),
    text = element_text(size = 10, face = "plain"),
    axis.title = element_text(size = 12, face = "plain", color = "black"),
    axis.text = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(colour = "black", angle = 0), 
    legend.title = element_text(size = 10, face = "plain", color = "black"),
    legend.text = element_text(size = 10, face = "plain", color = "black"),
    legend.background = element_blank(),
    legend.position = "right",
    axis.line.y = element_line(color = "black", linetype = "solid", linewidth = 0.2),
    axis.line.x = element_line(color = "black", linetype = "solid", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.2, fill = NA)
  )
  
  pcoa <- ggplot(points, aes(x = PCoA1, y = PCoA2, color = cg, shape = sg)) + 
    geom_point(size = point_size) +
    labs(x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
         y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = ""),
         title = paste(title, "R2 = ", round(adonis_res$R2[1], digits = 2), "P = ", adonis_res$`Pr(>F)`[1])) +
    theme_classic() +
    mytheme
  
  if (text_name == TRUE) {
    pcoa <- pcoa + geom_text(label = rownames(points), size = point_size, nudge_y = 0.02)
  }
  
  if (ellipse == TRUE) {
    pcoa <- pcoa + stat_ellipse(aes(fill = cg), geom = 'polygon', linetype = "solid", level = 0.6, alpha = 0.1, show.legend = TRUE)
  }
  
  return(list(bray, adonis_res, pcoa, tab2))
}
