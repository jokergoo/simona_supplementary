library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(cowplot)
library(simona)

wrap_plot = function(p, 
    margin_top = unit(0, "pt"),
    margin_left = unit(0, "pt"),
    margin_bottom = unit(0, "pt"), 
    margin_right = unit(0, "pt")) {
    if(inherits(p, "ggplot")) {
        p = grid.grabExpr(print(p))
    }

    vp = viewport(x = margin_left, y = unit(1, "npc") - margin_top, 
        width = unit(1, "npc") - margin_right - margin_left, 
        height = unit(1, "npc") - margin_bottom - margin_top,
        just = c("left", "top"))

    gTree(children = gList(p), vp = vp)
}
