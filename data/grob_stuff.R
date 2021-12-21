g_legend  = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

align.grobs = function(plots){
  grobs = list()
  widths = list()
  for (i in 1:length(plots)){
    grobs[[i]] = ggplotGrob(plots[[i]])
    widths[[i]] = grobs[[i]]$widths[2:5]
  }
  maxwidth = do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] = as.list(maxwidth)
  }
  grobs
}