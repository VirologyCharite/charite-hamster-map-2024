
rm(list = ls())

library(Racmacs)

### PRNT map
# Read the titer table
titers50 <- read.titerTable("data/220511/table-PRNT-220511.csv")

# Read in information about serum groups
sr_groups_info <- read.csv("data/metadata/sr_groups.csv", stringsAsFactors = FALSE)

# Make the map
map_full <- acmap(
  ag_names = rownames(titers50),
  sr_names = colnames(titers50),
  titer_table = titers50
)

dilutionStepsize(map_full) <- 0

# Set serum groups
sr_groups <- sr_groups_info$Serum.group[match(srNames(map_full), sr_groups_info$Serum)]
srGroups(map_full) <- factor(
  sr_groups,
  levels = c(
    "Pfizer",
    "AstraZeneca",
    "WT convalescent",
    "B.1.1.7 convalescent",
    "B.1.351 convalescent",
    "BA.1 convalescent"
  )
)

# Set antigen colors
ag_colors_info <- read.csv("data/metadata/ag_colors.csv", stringsAsFactors = FALSE)
ag_colors <- ag_colors_info$Color[match(agNames(map_full), ag_colors_info$Antigen)]
agFill(map_full) <- ag_colors

# Set serum colors
sr_colors_info <- read.csv("data/metadata/sr_colors.csv", stringsAsFactors = FALSE)
sr_colors <- sr_colors_info$Color[match(srGroups(map_full), sr_colors_info$Serum)]
srOutline(map_full) <- sr_colors

# Set styles
srOutlineWidth(map_full) <- 1
srSize(map_full) <- 10
agSize(map_full) <- 18

# Subset to only detectable titers and antigens (remove BA.1)
ndetectable <- colSums(titerTable(map_full) != "*" & !grepl("<", titerTable(map_full)) & !grepl(">", titerTable(map_full)))
map_subset <- subsetMap(map_full, sera = ndetectable > 2, antigens = c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2"))

# Optimize the maps
map_full <- optimizeMap(
  map_full,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = "none",
  options = list(ignore_disconnected = TRUE)
)
ptDrawingOrder(map_full) <- rev(ptDrawingOrder(map_full))

map_full <- rotateMap(map_full, -25)
save.acmap(map_full, 'data/220511/map-full-PRNT-220511.ace')
