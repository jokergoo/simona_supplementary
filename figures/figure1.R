

library(simona)
dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db")


lt = lapply(all_term_IC_methods(), function(method) {
    term_IC(dag, method)
})
names(lt) = all_term_IC_methods()
df = as.data.frame(lt)

cor = cor(df, use = "pairwise.complete.obs", method = "spearman")
hc = hclust(as.dist(1 - cor))

library(ComplexHeatmap)
p1 = grid.grabExpr(draw(Heatmap(cor, name = "correlation", cluster_rows = hc, cluster_columns = hc,
    column_title = "A) Correlation between ICs, on GO BP")))

group = c("IC_annotation" = 1, "IC_Seddiqui_2010" = 1, "IC_Zhang_2006" = 1,
          "IC_offspring" = 1, "IC_Seco_2004" = 1, "IC_universal" = 2,
          "IC_Wang_2007" = 2, "IC_Meng_2012" = 3, "IC_Sanchez_2011" = 3,
          "IC_height" = 3, "IC_Zhou_2008" = 3)


library(ggrepel)
library(ggplot2)
loc = cmdscale(as.dist(1 - cor))
loc = as.data.frame(loc)
colnames(loc) = c("x", "y")
loc$method = rownames(loc)

loc$group = group[rownames(loc)]

p2 = ggplot(loc, aes(x, y, label = method, col = factor(group))) + 
    geom_point() + 
    geom_text_repel(show.legend = FALSE, size = 3) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Group") +
    ggtitle("B) MDS based on the correlation between ICs")


#### sim
set.seed(123)
ic = term_IC(dag, method = "IC_annotation")
ic = ic[!is.na(ic)]
go_id = sample(names(ic), 500)


lt = lapply(all_term_sim_methods(), function(method) {
    term_sim(dag, go_id, method)
})
names(lt) = all_term_sim_methods()


df = as.data.frame(lapply(lt, function(x) x[lower.tri(x)]))

cor = cor(df, use = "pairwise.complete.obs")

ind = which(colnames(df) %in% c("Sim_Jiang_1997", 
    "Sim_Dice", "Sim_Kappa", "Sim_Jaccard", "Sim_Overlap", 
    "Sim_AIC_2014", "Sim_universal", "Sim_HRSS_2013"))

cor2 = cor[-ind, -ind]

p3 = grid.grabExpr(draw(Heatmap(cor2, name = "correlation", 
	column_title = "C) Correlation between semantic similarities, on GO BP")))

group = c("Sim_Pekar_2002" = 1, "Sim_Stojanovic_2001" = 1, "Sim_WP_1994" = 1,
          "Sim_Shenoy_2012" = 1, "Sim_Li_2003" = 1, "Sim_Wang_edge_2012" = 1,
          "Sim_SSDD_2013" = 2, "Sim_RSS_2013" = 2, "Sim_Zhong_2002" = 2,
          "Sim_Slimani_2006" = 2, "Sim_Shen_2010" = 3, "Sim_Zhang_2006" = 3,
          "Sim_EISI_2015" = 3, "Sim_XGraSM_2013" = 3, "Sim_Lin_1998" = 3,
          "Sim_Resnik_1999" = 3, "Sim_FaITH_2010" = 3, "Sim_Relevance_2006" = 3,
          "Sim_SimIC_2010" = 3, "Sim_Wang_2007" = 4, "Sim_Ancestor" = 4,
          "Sim_GOGO_2018" = 4, "Sim_AlMubaid_2006" = 4, "Sim_Rada_1989" = 4,
          "Sim_Leocock_1998" = 4, "Sim_Resnik_edge_2005" = 4)

loc = cmdscale(as.dist(1-cor2))
loc = as.data.frame(loc)
colnames(loc) = c("x", "y")
loc$method = rownames(loc)

loc$group = group[rownames(loc)]

p4 = ggplot(loc, aes(x, y, label = method, col = factor(group))) + 
    geom_point() + 
    geom_text_repel(show.legend = FALSE, size = 3) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Group") +
    ggtitle("D) MDS based on the correlation between similarities")



library(cowplot)

pdf("figure1.pdf", width = 12, height = 12)
print(plot_grid(p1, p2, p3, p4, nrow = 2))
dev.off()