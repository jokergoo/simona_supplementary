
source("lib.R")


library(ggplot2)

library(simona)
dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db")


set.seed(123)
ic = term_IC(dag, method = "IC_annotation")
ic = ic[!is.na(ic)]
go_id = sample(names(ic), 500)

term1 = MICA_term(dag, go_id)
term1 = term1[lower.tri(term1)]

term2 = LCA_term(dag, go_id)
term2 = term2[lower.tri(term2)]

l = term1 != term2

data = data.frame(x = 2, group = c("same", "different"), n = c(sum(!l), sum(l)))
data$label = paste0(round(data$n/sum(data$n)*100, 1), "%")
data$label[1] = paste0(data$label[1], "\nsame")
data$label[2] = paste0(data$label[2], "\ndifferent")
p1 = ggplot(data, aes(x = x, y = n, fill = group)) +
	geom_col() +
        geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
	coord_polar("y", start = 0) +
        xlim(c(1, 2.5)) + labs(y = "Number of term pairs") +
	ggtitle("A) Compare MICA and LCA") +
        theme(axis.line=element_blank(),
          # axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          # axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")


term1 = term1[l]
term2 = term2[l]

term_unique = unique(c(term1, term2))
dist = shortest_distances_via_NCA(dag, term_unique) # a symmetric matrix with terms as names
x = numeric(length(term1))
for(i in seq_along(term1)) {
    x[i] = dist[term1[i], term2[i]]
}

data = data.frame(x = x)
p2 = ggplot(data, aes(x = x)) + 
	geom_bar() +
	labs(x = "Undirected distance bewteem MICA and LCA", y = "Count") +
        ggtitle("B) On term pairs with different MICA and LCA")


mat1 = term_sim(dag, go_id, method = "Sim_Lin_1998")
mat2 = term_sim(dag, go_id, method = "Sim_WP_1994")

od = hclust(dist(mat1))$order

library(ComplexHeatmap)

ht_list = 
Heatmap(mat1, name = "similarity_Lin", 
        show_row_names = FALSE, show_column_names = FALSE,
        row_order = od, column_order = od, use_raster = TRUE, raster_quality = 2,
        column_title = "C) Sim_Lin_1998, based on MICA/IC_annotation") +
Heatmap(mat2, name = "similarity_WP", 
        show_row_names = FALSE, show_column_names = FALSE,
        row_order = od, column_order = od, use_raster = TRUE, raster_quality = 2,
        column_title = "D) Sim_WP_1994, based on LCA/depth")

p3 = grid.grabExpr(draw(ht_list, gap = unit(5, "mm")), width = 10, height = 5)

library(cowplot)
pdf("figure4.pdf", width = 10, height = 10)
print(plot_grid(plot_grid(
        p1, p2, nrow = 1), 
        wrap_plot(p3, margin_left = unit(20, "pt")), nrow = 2))
dev.off()
