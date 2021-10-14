# -- Libraries
library(scales)
library(tidyverse)
library(lubridate)
library(excessmort)
library(directlabels)

# -- Puerto Rico hurricane dates
hurricane_dates  <- c(Hugo    = make_date(1989, 9, 18),
                      Georges = make_date(1998, 9, 21),
                      Maria   = make_date(2017, 9, 20))


# Excluding  1) years (starting in 7/1) that include events, 2) 2020 and 3) some outliers in 2001 -------------------
exclude_dates <- c(make_date(1989, 7, 1) + 0:365,
                   make_date(1998, 7, 1) + 0:365,
                   make_date(2017, 7, 1) + 0:365,
                   make_date(2014, 7, 1) + 0:365,
                   seq(make_date(2001, 1, 1), make_date(2001, 1, 15), by = "day"),
                   seq(make_date(2020, 1, 1), today(), by = "day"))


# define control regions --------------------------------------------------
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")

# collapse age groups -----------------------------------------------------
data("puerto_rico_counts")
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)

# -- My color palette
my_palette        <- c("#456e9d", "#f08127", "#6baeaa", "#509745",  "#E69F00", "#a66f97", "#ff929d", "#D22B2B", "#252525", "#525252") # Old yellow: "#eac240"
names(my_palette) <- c("blue", "orange", "turquoise", "green", "yellow", "purple", "pink", "red", "black", "gray")

# -- Set up for figures
theme_set(theme_light(base_size   = 12, 
                      base_family = "Helvetica"))

# -- Modifying plot elements globally
theme_update(
  axis.ticks        = element_line(color = "grey92"),
  axis.ticks.length = unit(.5, "lines"),
  panel.grid.minor  = element_blank(),
  legend.title      = element_text(size = 12),
  legend.text       = element_text(color = "grey30"),
  legend.background = element_rect(color = "black", fill = "white"), #FBFCFC
  legend.key        = element_rect(fill = "white"),
  legend.direction  = "horizontal",
  legend.position   = "top",
  plot.title        = element_text(size = 18, face = "bold"),
  plot.subtitle     = element_text(size = 12, color = "grey30"),
  plot.caption      = element_text(size = 9, margin = margin(t = 15)),
  plot.background   = element_rect(fill="white", color = "white"), 
  panel.background  = element_rect(fill="white", color = "white"),
  strip.text        = element_text(face = "bold", color = "white"),
  strip.background  = element_rect(fill = "#252525"))

