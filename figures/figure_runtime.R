
source("lib.R")

library(ggplot2)
library(RColorBrewer)
library(GetoptLong)

###
load("../runtime/runtime_OBOFoundry_all.RData")

l = sapply(names(lt), function(nm) {
    max(lt[[nm]]$k)/df[nm, "n_terms"] > 0.8
})

lt = lt[l]

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


lines = rbind(data.frame(x = seq(0, 1, length = 2), y = seq(0, 1, length = 2), complexity = "O(n)"),
              data.frame(x = seq(0, 1, length = 100), y = seq(0, 1, length = 100)^2, complexity = "O(n^2)"),
              data.frame(x = seq(0, 1, length = 100), y = seq(0, 1, length = 100)^3, complexity = "O(n^3)"))
p1 = ggplot(df2, aes(x, y, col =)) +
    geom_line(aes(group = ontology), alpha = 0.25, show.legend = FALSE) +
    geom_line(data = lines, mapping = aes(x = x, y = y, col = complexity), lty = 2) +
    labs(x = "Numbers of terms, scaled", y = "runtime, scaled") +
    ggtitle("Compare runtime performance of OBOFoundry ontologies") +
    theme(legend.position = c(0.85, 0.16))


rel_diff = function(x, y) {
    if(missing(y)) {
        y = x[[2]]
        x = x[[1]]
    }

    x = x/max(x); x = c(0, x)
    y = y/max(y); y = c(0, y)

    fit = loess(y ~ x, span = 0.5)

    x2 = seq(0, 1, length = 1000)
    y2 = predict(fit, x2)

    area = sum(1/1000 * y2)

    0.5 - area
}


df$rel_diff = sapply(lt, rel_diff)

library(ggrepel)

p2 = ggplot(df, aes(x = n_terms, y = rel_diff)) +
    geom_point() +
    scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6), labels = c("1K", "10K", "100K", "1M")) +
    geom_hline(data = data.frame(yintercept = c(0, 1/6, 1/4), complexity = c("O(n)", "O(n^2)", "O(n^3)")), 
        mapping = aes(yintercept = yintercept, col = complexity), lty = 2) +
    lims(y = c(-0.3, 0.3)) +
    labs(x = "Numbers of terms", y = "Relative difference to linear time complexity") +
    ggtitle("Compare time complexity on OBOFoundry ontologies") +
    geom_text_repel(data = df[df$n_terms > 100000 & df$rel_diff < 0.1, ], 
        mapping = aes(x = n_terms, y = rel_diff, label = id), col = 1) +
    theme(legend.position = c(0.13, 0.16))


get_t_by_k = function(k) {
    sapply(lt, function(df) {
        i0 = which(df$k == k)
        if(length(i0)) {
            df$t[i0]
        } else {
            ind1 = which(df$k < k)
            if(length(ind1) == 0) { # k is smaller than all df$k
                return(NA)
            }
            i1 = max(ind1)
            ind2 = which(df$k > k)
            if(length(ind2) == 0) { # k is larger than all df$k
                return(NA)
            }
            i2 = min(ind2)

            x = c(0, df$k)
            y = c(0, df$t)
            fit = loess(y ~ x, span = 0.5)
            predict(fit, k)
         
        }
    })
}

df$t500 = get_t_by_k(500)

p3 = ggplot(df, aes(x = n_terms, y = t500, col = avg_parents)) +
    geom_point() +
    scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")))+
    scale_x_log10() + scale_y_log10() +
    labs(x = "Size of the ontology / n", y = "Runtime (sec) / t") +
    ggtitle("Randomly sample 500 terms from each ontology") +
    theme(legend.position = c(0.85, 0.22))

n = nrow(df)
fit = lm(log(df$t500) ~ log(df$n_terms))
su = summary(fit)
coef = coef(fit)
f = function(x) {
    exp(log(x)*coef[2] + coef[1])
}
p3 = p3 + geom_line(data = data.frame(x = range(df$n_terms), y = f(range(df$n_terms))), mapping = aes(x = x, y= y), col = "black")
coef = round(coef, 2)
p3 = p3 + annotate(geom = "text", x = 5e3, y = 44, label = qq("log(t) = @{coef[2]}*log(n) - @{-coef[1]}"))

compare_runtime_by_k = function(k) {
    tt = get_t_by_k(k)
    n = length(tt)
    fit = lm(log(tt) ~ log(df$n_terms))
    su = summary(fit)
    p = pf(su$fstatistic[1], su$fstatisti[2], su$fstatisti[3], lower.tail = FALSE)
    print(p)
    coef(fit)[2]
}

k = c(500, 600, 700, 800, 900, 1000, 2000, 3000, 
      4000, 5000, 6000, 7000, 8000, 9000, 10000)
coef = sapply(k, compare_runtime_by_k)


p4 = ggplot(data.frame(k = k, coef = coef), aes(x = k, y = coef)) +
    geom_point() +
    labs(x = "Numbers of random terms to fix", y = "Slope coef") +
    ggtitle("Linear regressions of runtime on the ontolgy size")


library(cowplot)

pdf("figure2.pdf", width = 10, height = 10)
print(plot_grid(
    p1, 
    p2, 
    wrap_plot(p3, margin_bottom = unit(5, "mm")), 
    wrap_plot(p4, margin_bottom = unit(5, "mm")), 
    nrow = 2
))
dev.off()
