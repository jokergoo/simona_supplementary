
source("lib.R")

library(GetoptLong)

library(cola)
data(golub_cola) 
res = golub_cola["ATC:skmeans"]

library(hu6800.db)
x = hu6800ENTREZID
mapped_probes = mappedkeys(x)
id_mapping = unlist(as.list(x[mapped_probes]))

lt = functional_enrichment(res, k = 3, id_mapping = id_mapping)

go_tb = lt[[1]]

sig_go_ids = go_tb$ID[go_tb$p.adjust < 0.01]
p.adjust = go_tb$p.adjust[go_tb$p.adjust < 0.01]

node_size_fun = function(x, range = c(2, 10)) {
	s = (range[2] - range[1])/(quantile(x, 0.95) - min(x)) * (x - min(x)) + range[1]
	s[s > range[2]] = range[2]
	s
}

lgd = Legend(title = "p.adjust", at = -log10(c(0.01, 0.001, 0.0001)), 
	labels = c("0.01", "0.001", "0.0001"), type = "points",
	size = unit(node_size_fun(-log10(c(0.01, 0.001, 0.0001))), "pt"))

dag = create_ontology_DAG_from_GO_db()


node_size = rep(2, dag_n_terms(dag))
names(node_size) = dag_all_terms(dag)
node_size[sig_go_ids] = node_size_fun(-log10(p.adjust))
pdf("figure1.pdf", width = 8, height = 5.2)
pushViewport(viewport(y = 0, height = unit(1, "npc") - unit(0.8, "cm"), just = "bottom"))
dag_circular_viz(dag, partition_by_size = round(dag_n_terms(dag)/5), 
	highlight = sig_go_ids,
	node_size = node_size,
	edge_transparency = 0.92, other_legends = lgd,
	legend_labels_max_width = 50,
	use_raster = TRUE, newpage = FALSE)
grid.text(qq("Gene Ontology, Biological Process, @{dag_n_terms(dag)} terms"), 
	y = unit(1, "npc") + unit(0.5, "mm"), just = "bottom", gp = gpar(fontsize = 14))
popViewport()
dev.off()
