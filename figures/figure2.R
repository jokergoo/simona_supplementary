


source("lib.R")

library(ggplot2)


library(simona)
dag = simona::create_ontology_DAG_from_GO_db(namespace = "BP", org_db = "org.Hs.eg.db", 
	relations = c("part_of", "regulates"))
ic1 = simona::term_IC(dag, method = "IC_annotation")

go_id1 = names(ic1[!is.na(ic1)])


library(ontologyIndex)
onto = ontologyIndex::get_ontology("../compare_tools/go-basic.obo", 
	propagate_relationships = c("is_a", "part_of", "regulates", "negatively_regulates", "positively_regulates"))

library(org.Hs.eg.db)
tb = toTable(org.Hs.egGO)
tb = tb[tb$Ontology == "BP", ]
# this ensures each GO term has genes annotated
lt = split(tb$go_id, tb$gene_id)

ic2 = ontologyIndex::get_term_info_content(onto, lt)
go_id2 = names(ic2)


library(GOSemSim)
data = GOSemSim::godata("org.Hs.eg.db", ont = "BP")
ic3 = data@IC[is.finite(data@IC)]
go_id3 = names(ic3)



library(GOSim)
GOSim::setOntology("BP")
GOSim::calcICs()
load("ICsBPhumanall.rda") # it loads the object `IC`
ic4 = IC[is.finite(IC)]
go_id4 = names(ic4)



go_id_common = intersect(intersect(go_id1, go_id2), intersect(go_id3, go_id4))
df = data.frame(
	simona        = ic1[go_id_common],
	ontologyIndex = ic2[go_id_common],
	GOSemSim      = ic3[go_id_common],
	GOSim         = ic4[go_id_common]
)


p1 = ggplot(df, aes(x = simona, y = ontologyIndex)) +
	geom_point(size = 1) +
	ggtitle("A) Compare ICs between simona / ontologyIndex")
p2 = ggplot(df, aes(x = simona, y = GOSemSim)) +
	geom_point(size = 1) + geom_abline(slope= 1, intercept=0, col = 2, lty = 2) +
	scale_x_continuous(limits = c(0, 12.5), breaks = seq(0, 12, by = 4)) +
	scale_y_continuous(limits = c(0, 12.5), breaks = seq(0, 12, by = 4)) +
	ggtitle("B) Compare ICs between simona / GOSemSim")
p3 = ggplot(df, aes(x = GOSim, y = GOSemSim)) +
	geom_point(size = 1) +
	scale_x_continuous(limits = c(0, 12.5), breaks = seq(0, 12, by = 4)) +
	scale_y_continuous(limits = c(0, 12.5), breaks = seq(0, 12, by = 4)) +
	ggtitle("C) Compare ICs between GOSim / GOSemSim")

#########################

library(ComplexHeatmap)

load("../compare_tools/compare_sim_mat.RData")

library(circlize)
# set the same color mapping to the four heatmaps
col_fun = colorRamp2(c(0, 0.3, 0.6), c("blue", "white", "red"))

od = hclust(dist(m1))$order


p4 = grid.grabExpr(
	draw(Heatmap(m1, name = "similarity", col = col_fun,
		show_row_names = FALSE, show_column_names = FALSE,
		row_order = od, column_order = od, show_heatmap_legend = FALSE, 
		use_raster = TRUE, raster_quality = 2,
		column_title = "D) Similaritysimona, random 500 GO terms")
	), width = 4, height = 4
)
p5 = grid.grabExpr(
	draw(Heatmap(m2, name = "similarity2", col = col_fun,
		show_row_names = FALSE, show_column_names = FALSE,
		row_order = od, column_order = od, show_heatmap_legend = FALSE, 
		use_raster = TRUE,raster_quality = 2,
		column_title = "E) ontologySimilarity, random 500 GO terms")
	), width = 4, height = 4
)
p6 = grid.grabExpr(
	draw(Heatmap(m3, name = "similarity3", col = col_fun, 
		show_row_names = FALSE, show_column_names = FALSE,
		row_order = od, column_order = od, show_heatmap_legend = FALSE, 
		use_raster = TRUE,raster_quality = 2,
		column_title = "F) GOSemSim, random 500 GO terms")
	), width = 4, height = 4
)
p7 = grid.grabExpr(
	draw(Heatmap(m4, name = "Similarity", col = col_fun, 
		show_row_names = FALSE, show_column_names = FALSE, 
		row_order = od, column_order = od,
		use_raster = TRUE, raster_quality = 2,
		column_title = "G) GOSim, random 500 GO terms")
	), width = 4, height = 4
)

####

load("../compare_tools/compare_tools.RData")

df1 = rbind(data.frame(k = k, runtime = t1, tool = "simona"),
	       data.frame(k = k, runtime = t2, tool = "ontologySimilarity"),
	       data.frame(k = k, runtime = t3, tool = "GOSemSim"), 
	       data.frame(k = k, runtime = t4, tool = "GOSim"))
df1$tool = factor(df1$tool, levels = c("simona", "ontologySimilarity", "GOSemSim", "GOSim"))
df1 = df1[!is.na(df1$runtime), ]


p8 = ggplot(df1, aes(x = k, y = runtime, col = tool)) +
	geom_line()	+ 
	labs(x = "Numbers of terms", y = "Runtime (sec)", col = "Tool") +
	ggtitle("H) Compare runtime performance")

p9 = grob()


library(cowplot)

pdf("figure2.pdf", width = 13, height = 13)
print(plot_grid(
	plot_grid(p1, p2, p3, nrow = 1),
	plot_grid(
		wrap_plot(p4, margin_left = unit(20, "pt")), 
		wrap_plot(p5, margin_left = unit(20, "pt")), 
		wrap_plot(p6, margin_left = unit(20, "pt")), 
		nrow = 1),
	plot_grid(
		wrap_plot(p7, margin_left = unit(20, "pt"), margin_right = unit(-50, "pt")), 
		wrap_plot(p8, margin_left = unit(60, "pt"), margin_right = unit(-170, "pt"), margin_bottom = unit(-20, "pt")), 
		p9, nrow = 1),
	p9,
	nrow = 4, rel_heights = c(1, 1, 1, 0.05)))
dev.off()

