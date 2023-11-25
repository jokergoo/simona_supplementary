
source("lib.R")

library(simona)
dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db")


lt = lapply(all_term_IC_methods(), function(method) {
    term_IC(dag, method)
})
names(lt) = all_term_IC_methods()
df = as.data.frame(lt)

cor_ic = cor(df, use = "pairwise.complete.obs")
hc_ic = hclust(as.dist(1 - cor_ic))

group_ic = c("IC_height" = "1",
          "IC_Sanchez_2011" = "1",
          "IC_Zhou_2008" = "1",
          "IC_Meng_2012" = "1",
          "IC_Seddiqui_2010" = "2",
          "IC_Seco_2004" = "2",
          "IC_offspring" = "2",
          "IC_Zhang_2006" = "2",
          "IC_annotation" = "2",
          "IC_Wang_2007" = "others",
          "IC_universal" = "others")


library(ComplexHeatmap)
p1 = grid.grabExpr(
    draw(Heatmap(cor_ic, name = "correlation", cluster_rows = hc_ic, cluster_columns = hc_ic,
        left_annotation = rowAnnotation(group = group_ic[rownames(cor_ic)], 
            col = list(group = c("1" = "#F8766D", "2" = "#00BA38", "others" = "#619CFF"))),
        top_annotation = HeatmapAnnotation(group = group_ic[rownames(cor_ic)],show_legend = FALSE,
            col = list(group = c("1" = "#F8766D", "2" = "#00BA38", "others" = "#619CFF"))),
        row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
        show_row_dend = FALSE, show_column_dend = FALSE,
        row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        column_title = "A) IC correlation, GO BP"))
)



library(ggrepel)
library(ggplot2)
loc = cmdscale(as.dist(1 - cor_ic))
loc = as.data.frame(loc)
colnames(loc) = c("x", "y")
loc$method = rownames(loc)

loc$group = group_ic[rownames(loc)]

p2 = ggplot(loc, aes(x, y, label = method, col = factor(group))) + 
    geom_point() + 
    geom_text_repel(show.legend = FALSE, size = 3) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Group") +
    ggtitle("B) MDS based on IC correlations")


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

cor_sim = cor(df, use = "pairwise.complete.obs")

ind = which(colnames(df) %in% c("Sim_Jiang_1997", "Sim_HRSS_2013", "Sim_universal",
    "Sim_Dice", "Sim_Kappa", "Sim_Jaccard", "Sim_Overlap"))
cor_sim2 = cor_sim[-ind, -ind]
df2 = df[, -ind]



group_sim = c("Sim_Shen_2010" = 1,
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


hc_sim = hclust(as.dist(1 - cor_sim2))
p3 = grid.grabExpr(
    draw(Heatmap(cor_sim2, name = "correlation", cluster_rows = hc_sim, cluster_columns = hc_sim,
        left_annotation = rowAnnotation(group = group_sim[rownames(cor_sim2)], 
            col = list(group = c("1" = "#F8766D", "2" = "#A3A500", "3" = "#00BF7D", "4" = "#00B0F6", "5" = "#E76BF3"))),
        top_annotation = HeatmapAnnotation(group = group_sim[rownames(cor_sim2)], show_legend = FALSE,
            col = list(group = c("1" = "#F8766D", "2" = "#A3A500", "3" = "#00BF7D", "4" = "#00B0F6", "5" = "#E76BF3"))),
        row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
        show_row_dend = FALSE, show_column_dend = FALSE,
        row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        column_title = "C) Semantic similarity correlation, GO BP"))
)


loc = cmdscale(as.dist(1-cor_sim2))
loc = as.data.frame(loc)
colnames(loc) = c("x", "y")
loc$method = rownames(loc)

loc$group = group_sim[rownames(loc)]

p4 = ggplot(loc, aes(x, y, label = method, col = factor(group))) + 
    geom_point() + 
    geom_text_repel(show.legend = FALSE, size = 3) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Group") +
    ggtitle("D) MDS based on similarity correlations")

p4 = grid.grabExpr(
    print(p4), width = 6, height = 5.5
)

library(cowplot)

pdf("figure3.pdf", width = 12, height = 11)
print(plot_grid(
    wrap_plot(p1, margin_right = unit(50, "pt"), margin_bottom = unit(17, "pt")), 
    wrap_plot(p2, margin_bottom = unit(20, "pt")), 
    p3, 
    wrap_plot(p4, margin_bottom = unit(1.3, "cm"), margin_right = unit(15, "pt")), 
    nrow = 2, rel_widths = c(1.1, 1), rel_heights = c(0.9, 1)
))
dev.off()
