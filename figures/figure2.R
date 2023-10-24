

library(ComplexHeatmap)

load("compare_sim_mat.RData")

library(circlize)
# set the same color mapping to the four heatmaps
col_fun = colorRamp2(c(0, 0.3, 0.6), c("blue", "white", "red"))

od = hclust(dist(m1))$order

ht_list = 
Heatmap(m1, name = "similarity", col = col_fun,
		show_row_names = FALSE, show_column_names = FALSE,
		row_order = od, column_order = od, use_raster = TRUE, raster_quality = 4,
		column_title = "simona, random 500 GO terms") +
Heatmap(m2, name = "similarity2", col = col_fun,
		show_row_names = FALSE, show_column_names = FALSE,
		row_order = od, column_order = od, show_heatmap_legend = FALSE, use_raster = TRUE,raster_quality = 4,
		column_title = "ontologySimilarity, random 500 GO terms") +
Heatmap(m3, name = "similarity3", col = col_fun, 
		show_row_names = FALSE, show_column_names = FALSE,
		row_order = od, column_order = od, show_heatmap_legend = FALSE, use_raster = TRUE,raster_quality = 4,
		column_title = "GOSemSim, random 500 GO terms") +
Heatmap(m4, name = "similarity4", col = col_fun, 
		show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE, use_raster = TRUE,raster_quality = 4,
		row_order = od, column_order = od,
		column_title = "GOSim, random 500 GO terms")

p1 = grid.grabExpr(draw(ht_list), width = 10, height = 10)


#########################

load("../compare_tools/compare_tools.RData")

df1 = rbind(data.frame(k = k, runtime = t1, tool = "simona"),
	       data.frame(k = k, runtime = t2, tool = "ontologySimilarity"),
	       data.frame(k = k, runtime = t3, tool = "GOSemSim"), 
	       data.frame(k = k, runtime = t4, tool = "GOSim"))
df1$tool = factor(df1$tool, levels = c("simona", "ontologySimilarity", "GOSemSim", "GOSim"))
df1 = df1[!is.na(df1$runtime), ]
df1 = df1[df1$runtime < 400, ]

library(ggplot2)

p2 = ggplot(df1, aes(x = k, y = runtime, col = tool)) +
	geom_line()	+
	labs(x = "Numbers of terms", y = "runtime (sec)", col = "Tool") +
	ggtitle("Compare runtime of four tools") +
	theme(legend.position = c(0.4, 0.99), legend.justification = c(0.5, 1))

load("../runtime/runtime_OBOFoundry_all.RData")

for(nm in names(lt)) {
	x = lt[[nm]]$k
    y = lt[[nm]]$t
    x = x/max(x)
    y = y/max(y)

    lt[[nm]]$x = x
    lt[[nm]]$y = y
	lt[[nm]]$ontology = nm
}
df2 = do.call(rbind, lt)

p3 = ggplot(df2, aes(x, y, col = ontology)) +
	geom_line(show.legend = FALSE) +
	scale_color_manual(values = rep("#00000040", 99)) +
	geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
	labs(x = "Numbers of terms, scaled", y = "runtime, scaled") +
	ggtitle("Compare runtime performance of OBOFoundry ontologies")


load("../runtime/runtime_OBOFoundry_all.RData")

rel_diff = function(x, y) {
    if(missing(y)) {
        y = x[[2]]
        x = x[[1]]
    }

    od = order(x)
    x = x[od]
    y = y[od]
    n = length(x)

    x = x/max(x)
    y = y/max(y)

    area = sum( 0.5*(x[2:n] - x[2:n - 1])*(y[2:n] + y[2:n - 1]) )

    0.5 - area
}


df$rel_diff = sapply(lt, rel_diff)

library(ggrepel)

p4 = ggplot(df, aes(x = n_terms, y = rel_diff)) +
	geom_point() +
	scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6), labels = c("1K", "10K", "100K", "1M")) +
	geom_hline(yintercept = 0, col = 4, lty = 2) + 
	labs(x = "Numbers of terms", y = "Relative difference to linear time complexity") +
	ggtitle("Compare time complexity on OBOFoundry ontologies") +
	geom_text_repel(data = df[df$n_terms > 100000 & df$rel_diff < 0.1, ], 
		mapping = aes(x = n_terms, y = rel_diff, label = id), col = 2)


library(cowplot)

pdf("figure2.pdf", width = 12*1.2, height = 7*1.2)
print(plot_grid(p1, plot_grid(p2, p3, p4, nrow = 1), nrow = 2, rel_heights = 3.4))
dev.off()

