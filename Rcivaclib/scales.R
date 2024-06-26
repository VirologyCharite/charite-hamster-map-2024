log2foldchange <- function(x) {

  xabs <- abs(x)
  foldchange <- 2^xabs
  foldchange[x < 0] <- -foldchange[x < 0]
  as.character(round(foldchange, 1))

}


titerplot_theme <- function(){

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "grey40"
    ),
    panel.grid.major = element_line(
      linetype = "solid",
      colour = rgb(0,0,0,0.05)
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      # face = "bold",
      size = 10
    )
  )

}


agFillScale <- function(map) {
  ag_fill <- agFill(map)
  names(ag_fill) <- agNames(map)
  ag_fill
}
