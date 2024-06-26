
rm(list = ls())

library(Racmacs)
library(tidyverse)

# Read in colour information
ag_colors_info <- read.csv("data/metadata/ag_colors.csv", stringsAsFactors = FALSE)
ag_names_info <- read.csv("data/metadata/ag_names.csv", stringsAsFactors = FALSE)
sr_groups_info <- read.csv("data/metadata/sr_groups.csv", stringsAsFactors = FALSE)
sr_colors_info <- read.csv("data/metadata/sr_colors.csv", stringsAsFactors = FALSE)

orig <- read.acmap('../../data/230420-hamster-retitrations/map-continuous-fixbottom-90-corrected.ace')

make_map <- function(titer_table, dilution_stepsize = 1) {
  discrete <- acmap(titer_table = titer_table)
  dilutionStepsize(discrete) <- dilution_stepsize
  discrete <- optimizeMap(discrete, 2, 500)
  agSize(discrete) <- 18
  srSize(discrete) <- 8
  srOutlineWidth(discrete) <- 2
  ptDrawingOrder(discrete) <- rev(ptDrawingOrder(discrete))
  
  # Colour the antigens
  ag_colors <- ag_colors_info$Color[match(agNames(discrete), ag_colors_info$Antigen)]
  agFill(discrete) <- ag_colors
  
  # Set the antigen names
  agNames(discrete) <- ag_names_info$Name[match(agNames(discrete), ag_names_info$Antigen)]
  
  # Set the serum groups
  serum_groups <- sr_groups_info$Serum.group[match(srNames(discrete), sr_groups_info$Serum)]
  srGroups(discrete) <- factor(serum_groups, levels = c("D614G sera", "Alpha sera", "Delta sera",
                                                        "Beta sera", "B.1+E484K sera", "BA.1 sera",
                                                        "BA.2 sera", "BA.2-12 sera", "BA.5 sera", 
                                                        "XBB.2 sera"))
  sr_colors <- sr_colors_info$Color[match(srGroups(discrete), sr_colors_info$Serum)]
  srOutline(discrete) <- sr_colors
  
  discrete <- realignMap(discrete, orig)
  discrete <- rotateMap(discrete, 20)
  discrete
}

## Make the discrete maps
# For different sensitivity levels
titer_table <- read.titerTable('data/240123-hamster/titers-discrete-50.csv')
disc50 <- make_map(titer_table = titer_table)
save.acmap(disc50, 'data/240123-hamster/map-discrete-50.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-discrete-75.csv')
disc75 <- make_map(titer_table = titer_table)
save.acmap(disc75, 'data/240123-hamster/map-discrete-75.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-discrete-90.csv')
disc90 <- make_map(titer_table = titer_table)
save.acmap(disc90, 'data/240123-hamster/map-discrete-90.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-discrete-99.csv')
disc99 <- make_map(titer_table = titer_table)
save.acmap(disc99, 'data/240123-hamster/map-discrete-99.ace')

## Make continuous maps
# For different sensitivity levels
titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixtop-fixbottom-50.csv')
cont50 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont50, 'data/240123-hamster/map-continuous-fixtop-fixbottom-50.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixtop-fixbottom-75.csv')
cont75 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont75, 'data/240123-hamster/map-continuous-fixtop-fixbottom-75.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixtop-fixbottom-90.csv')
cont90 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont90, 'data/240123-hamster/map-continuous-fixtop-fixbottom-90.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixtop-fixbottom-99.csv')
cont99 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont99, 'data/240123-hamster/map-continuous-fixtop-fixbottom-99.ace')

# For different methods
titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixtop-90.csv')
cont90 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont90, 'data/240123-hamster/map-continuous-fixtop-90.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixbottom-90.csv')
cont90 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont90, 'data/240123-hamster/map-continuous-fixbottom-90.ace')

titer_table <- read.titerTable('data/240123-hamster/titers-continuous-90.csv')
cont90 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont90, 'data/240123-hamster/map-continuous-90.ace')

# PRNT90 fixbottom, adapted
titer_table <- read.titerTable('data/240123-hamster/titers-continuous-fixbottom-90-corrected.csv')
cont90 <- make_map(titer_table = titer_table, dilution_stepsize = 0)
save.acmap(cont90, 'data/240123-hamster/map-continuous-fixbottom-90-corrected.ace')

