# Alternative to the Venn diagram
library(UpSetR)
library(readxl)

setwd("I:/")

df <- read_excel("upSet.xlsx", col_names = TRUE)

listInput <- as.list(as.data.frame(df))

upset(fromList(listInput), nsets = 11, point.size = 3.5, line.size = 1.2,
      mainbar.y.label = "Pathway Intersections", sets.x.label = "Number of Pathways",
      text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 1.1), order.by = c("freq", "degree"),
      queries = list(list(query = intersects, params = list("miR1", "miR9"), color = "red", active = T)),
      empty.intersections = "on")

upset(fromList(listInput), nsets = 11, nintersects = 100, point.size = 3.5, line.size = 1.2,
      mainbar.y.label = "Pathway Intersections", sets.x.label = "Number of Pathways", group.by = "sets", cutoff = 7,
      text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 1.1))

intersect(as.vector(na.omit(df$miR1)), as.vector(na.omit(df$miR9)))
