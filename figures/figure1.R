

library(simona)
dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db")


lt = lapply(all_term_IC_methods(), function(method) {
    term_IC(dag, method)
})
names(lt) = all_term_IC_methods()
df = as.data.frame(lt)

cor = cor(df, use = "pairwise.complete.obs")
hc = hclust(as.dist(1 - cor))

library(ComplexHeatmap)
p1 = grid.grabExpr(
    draw(Heatmap(cor, name = "correlation", cluster_rows = hc, cluster_columns = hc,
        row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        column_title = "IC correlation, GO BP"))
)

group = c("IC_height" = "1",
          "IC_Sanchez_2011" = "1",
          "IC_Zhou_2008" = "1",
          "IC_Seddiqui_2010" = "2",
          "IC_Seco_2004" = "2",
          "IC_offspring" = "2",
          "IC_Zhang_2006" = "2",
          "IC_annotation" = "2",
          "IC_Meng_2012" = "others",
          "IC_Wang_2007" = "others",
          "IC_universal" = "others")


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
    ggtitle("MDS based on the correlation between ICs from different IC methods")


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

ind = which(colnames(df) %in% c("Sim_Jiang_1997", "Sim_HRSS_2013", "Sim_universal",
    "Sim_Dice", "Sim_Kappa", "Sim_Jaccard", "Sim_Overlap"))
cor2 = cor[-ind, -ind]
df2 = df[, -ind]

hc = hclust(as.dist(1 - cor2))
p3 = grid.grabExpr(
    draw(Heatmap(cor2, name = "correlation", cluster_rows = hc, cluster_columns = hc,
        row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        column_title = "Semantic similarity correlation, GO BP"))
)

group = c("Sim_Shen_2010" = 1,
          "Sim_Zhang_2006" = 1,
          "Sim_EISI_2015" = 1,
          "Sim_XGraSM_2013" = 1,
          "Sim_Resnik_1999" = 1,
          "Sim_Lin_1998" = 1,
          "Sim_FaITH_2010" = 1,
          "Sim_Relevance_2006" = 1,
          "Sim_SimIC_2010" = 1,
          "Sim_SSDD_2013" = 2,
          "Sim_RSS_2013" = 2,
          "Sim_Zhong_2002" = 2,
          "Sim_Slimani_2006" = 2,
          "Sim_Pekar_2002" = 3,
          "Sim_WP_1994" = 3,
          "Sim_Shenoy_2012" = 3,
          "Sim_Stojanovic_2001" = 3,
          "Sim_Li_2003" = 3,
          "Sim_Wang_edge_2012" = 3,
          "Sim_Wang_2007" = 4,
          "Sim_Ancestor" = 4,
          "Sim_AIC_2014" = 4,
          "Sim_GOGO_2018" = 4,
          "Sim_AlMubaid_2006" = 5,
          "Sim_Leocock_1998" = 5,
          "Sim_Rada_1989" = 5,
          "Sim_Resnik_edge_2005" = 5)


loc = cmdscale(as.dist(1-cor2))
loc = as.data.frame(loc)
colnames(loc) = c("x", "y")
loc$method = rownames(loc)

loc$group = group[rownames(loc)]

p4 = ggplot(loc, aes(x, y, label = method, col = factor(group))) + 
    geom_point() + 
    geom_text_repel(show.legend = FALSE, size = 3) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Group") +
    ggtitle("MDS based on the correlation between similarities from different term sim methods")



library(cowplot)

pdf("figure1.pdf", width = 12, height = 12)
print(plot_grid(p1, p2, p3, p4, nrow = 2))
dev.off()
