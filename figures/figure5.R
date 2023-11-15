
source("lib.R")


dag1 = import_owl("~/workspace/ontology/OBOFoundry/ncbitaxon/ncbitaxon.owl",
    remove_rings = TRUE, remove_cyclic_paths = TRUE)

dag2 = import_obo("~/workspace/ontology/OBOFoundry/ncit/ncit.obo",
    remove_rings = TRUE, remove_cyclic_paths = TRUE)

dag3 = import_owl("~/workspace/ontology/OBOFoundry/ogg/ogg.owl",
    remove_rings = TRUE, remove_cyclic_paths = TRUE)

dag4 = create_ontology_DAG_from_GO_db()


pdf("figure5.pdf", width = 20, height = 10)

fraction = 0.45
set.seed(123)
grid.newpage()
pushViewport(viewport(x = 0, y = 0.5, width = fraction, height = unit(0.5, "npc") - unit(5, "mm"), just = c("left", "bottom")))
dag_circular_viz(dag1, partition_by_size = round(dag_n_terms(dag1)/5), 
	node_transparency = 0.75, edge_transparency = 0.85, legend_labels_max_width = 25,
	use_raster = TRUE, newpage = FALSE)
grid.text(qq("NCBI organismal classification / ncbitaxon, @{dag_n_terms(dag1)} terms"), y = unit(1, "npc"), just = "top", gp = gpar(fontsize = 14))
popViewport()

pushViewport(viewport(x = fraction, y = 0.5, width = 1 - fraction, height = unit(0.5, "npc") - unit(5, "mm"), just = c("left", "bottom")))
dag_circular_viz(dag2, partition_by_size = round(dag_n_terms(dag2)/5), 
	node_transparency = 0.75, legend_labels_max_width = 50,
	use_raster = TRUE, newpage = FALSE)
grid.text(qq("NCI Thesaurus / ncit, @{dag_n_terms(dag2)} terms"), y = unit(1, "npc"), just = "top", gp = gpar(fontsize = 14))
popViewport()

pushViewport(viewport(x = 0, y = 0, width = fraction, height = unit(0.5, "npc") - unit(5, "mm"), just = c("left", "bottom")))
dag_circular_viz(dag3, partition_by_size = round(dag_n_terms(dag3)/5), 
	legend_labels_max_width = 25, edge_transparency = 0.75,
	use_raster = TRUE, newpage = FALSE)
grid.text(qq("The Ontology of Genes and Genomes / ogg, @{dag_n_terms(dag3)} terms"), y = unit(1, "npc"), just = "top", gp = gpar(fontsize = 14))
popViewport()

sig_go_ids = readRDS(system.file("extdata", "sig_go_ids.rds", package = "simona"))
pushViewport(viewport(x = fraction, y = 0, width = 1 - fraction, height = unit(0.5, "npc") - unit(5, "mm"), just = c("left", "bottom")))
dag_circular_viz(dag4, partition_by_size = round(dag_n_terms(dag4)/5), 
	highlight = sig_go_ids,
	legend_labels_max_width = 50,
	use_raster = TRUE, newpage = FALSE)
grid.text(qq("Gene Ontology, Biological process, @{dag_n_terms(dag4)} terms"), y = unit(1, "npc"), just = "top", gp = gpar(fontsize = 14))
popViewport()

dev.off()
