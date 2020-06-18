#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Preliminary assessment of the regional conservation status of snubfin
#'         dolphins (Orcaella heinsohni) in Western Australia
#'         
#'         Bouchet et al. (2020). Frontiers in Marine Science
#' 
#'                  --- 04: FIGURES & TABLES ---
#'                  
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------

# Set ggplot2 options

gg.opts <- theme(text = element_text(family = "Arial"),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 plot.title = element_text(size = 13, face = "bold"),
                 legend.title = element_text(size = 12), 
                 legend.text = element_text(size = 11, colour = "black"),
                 axis.text = element_text(size = 12, colour = "black"),
                 panel.grid.major = element_line(colour = 'lightgrey', size = 0.1),
                 axis.title.x = element_text(size = 13, margin = margin(t = 15, r = 0, b = 0, l = 0)),
                 axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 15, b = 0, l = 0)),
                 legend.position = "right")

# Tables ------------------------------------------------------------------

# | Table 1 ----------

# Number of sightings per dataset

snub %>% dplyr::group_by(dataset_ID) %>% dplyr::summarise(N = n())

# Effort (length of surveyed tracklines)

effort.km <- purrr::map(.x = tracks, .f = ~.x %>% sf::st_as_sf(.) %>% sf::st_length(.) %>% sum(.)/1000) %>% 
  reshape2::melt(.) %>% 
  dplyr::mutate(value = as.numeric(value)) %>% 
  dplyr::rename(track_length = value, dataset_ID = L1)


# Figures -----------------------------------------------------------------

load("data/gmap/gmap.RData")

# | Figure 1 ----------

map.kimb <- map_sightings(baselayer = gmap.kimb)
map.roebuck <- map_sightings(baselayer = gmap.roebuck)
map.cygnet <- map_sightings(baselayer = gmap.cygnet)
map.prince <- map_sightings(baselayer = gmap.prince)
map.cambridge <- map_sightings(baselayer = gmap.cambridge)

ggsave(filename = "out/Fig1_Kimb.pdf", plot = map.kimb, device = "pdf", dpi = 650, width = 8, height = 10)
ggsave(filename = "out/Fig1_Roebuck.pdf", plot = map.roebuck, device = "pdf", dpi = 650, width = 8, height = 10)
ggsave(filename = "out/Fig1_Cygnet.pdf", plot = map.cygnet, device = "pdf", dpi = 650, width = 8, height = 10)
ggsave(filename = "out/Fig1_Prince.pdf", plot = map.prince, device = "pdf", dpi = 650, width = 8, height = 10)
ggsave(filename = "out/Fig1_Cambridge.pdf", plot = map.cambridge, device = "pdf", dpi = 650, width = 8, height = 10)


#  | Figure 3 ----------

# Define custom colour ramp

color.ramp <- colorspace::sequential_hcl(n = 100, 
                                         h = c(285, 75), # hue
                                         c = c(30, NA, 73),
                                         l = c(20, 88), # luminosity
                                         power = c(2, 1.15), rev = FALSE)

# Bivariate kernel density plot 

fig3A <- raster::as.data.frame(bivariate.r, xy = TRUE) %>% 
  ggplot(data = .) + 
  geom_raster(aes(x = x, y = y, fill = layer)) + 
  xlab("Distance to river (km)") + 
  ylab("Depth (m)") + 
  scale_x_continuous(breaks = seq(0, 80, 10), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-100, 0, 10), expand = c(0,0)) +
  scale_fill_gradientn(name = "", breaks = seq(0, 1, 0.1), colors = pals::parula(100), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  gg.opts +
  theme(axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5))) +
  theme(legend.key.height = unit(1.45, "cm"), 
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", size = 0.8, fill = "transparent"))

# Boxplot

fig3B <- snub.dat %>% dplyr::filter(sighting_class=="secondary") %>% 
  droplevels() %>% 
  dplyr::mutate(dataset_ID = gsub(pattern = "OH_", replacement = "", x = dataset_ID)) %>% 
  ggplot(., aes(x = dataset_ID, y = incl.prob)) + 
  stat_boxplot(geom = 'errorbar', width = 0.25) +
  geom_boxplot(fill = "lightgrey", varwidth = TRUE) + gg.opts +
  xlab("Dataset ID") + 
  ylab("Inclusion probability") + 
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)))

# Combine both sub-figures

fig3 <- (fig3A | fig3B) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 18))

# # Save to file
# 
# ggsave(filename = "/Users/philippebouchet/Google Drive/Documents/projects/orcaella/AMMC/Papers/Distribution/FMS/Figures/Bouchet_etal_Fig3.png", plot = fig3, height = 4, width = 10)


#  | Figure 4 ----------

eoo.ci <- bci(eoo.mcp) # Calculate CI for MCP-EOO
aoo.ci <- bci(aoo) # Calculate CI for AOO

fig4A <- data.frame(eoo = eoo.mcp) %>% 
  ggplot(data = ., aes(x = eoo, y = 1000*..density..)) +
  geom_histogram(binwidth = 750, color = "white", fill = "grey75", size = 0.4, alpha = 1) +
  gg.opts + geom_density(col = "grey30", fill = "white", alpha = 0.4) +
  scale_x_continuous(label = comma, expand = c(0,0)) +
  scale_y_continuous(expand = c(0.03,0), breaks = seq(0, 0.15, 0.025)) +
  theme(axis.ticks = element_blank()) +
  bbc_style() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.ticks.x = ggplot2::element_line(color = "#cbcbcb"),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10))) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.2, 
               data = data.frame(x1 = eoo.ci[2], x2 = eoo.ci[3], y1 = 0.0, y2 = 0.0)) +
  geom_point(aes(x = x1, y = y1), colour = "black", size = 4,
               data = data.frame(x1 = eoo.ci[1], y1 = 0.0)) +
  ylab(expression(Density~("x"~10^3))) + 
  xlab("Extent of Occurrence (EOO)")

fig4B <- data.frame(aoo = aoo$all) %>% 
  ggplot(data = ., aes(x = aoo, y = 1000*..density..)) +
  geom_histogram(binwidth = 5, color = "white", fill = "grey75", size = 0.4, alpha = 1) +
  gg.opts + geom_density(col = "grey30", fill = "white", alpha = 0.4) +
  scale_x_continuous(label = comma, expand = c(0,0), breaks = seq(640, 760, 20)) +
  scale_y_continuous(expand = c(0.03,0), breaks = seq(0,30,5)) +
  theme(axis.ticks = element_blank()) +
  bbc_style() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.ticks.x = ggplot2::element_line(color = "#cbcbcb"),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.2, 
               data = data.frame(x1 = aoo.ci[2], x2 = aoo.ci[3], y1 = 0, y2 = 0)) +
  geom_point(aes(x = x1, y = y1), colour = "black", size = 4,
             data = data.frame(x1 = aoo.ci[1], y1 = 0.0)) +
  ylab(expression(Density~("x"~10^3))) + 
  xlab("Area of occupancy (AOO)")

fig4 <- (fig4A | fig4B) + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 18))

# # Save to file
# 
# ggsave(filename = "/Users/philippebouchet/Google Drive/Documents/projects/orcaella/AMMC/Papers/Distribution/FMS/Figures/Bouchet_etal_Fig4.png", plot = fig4, height = 4.5, width = 10)

# Supplementary -----------------------------------------------------------------

#  | Appendix 1 ----------

eoo.alpha.ci <- bci(eoo.alpha) # Calculate CI for EOO

fig.A1 <- data.frame(eoo = eoo$all) %>% 
  ggplot(data = ., aes(x = eoo, y = 1000*..density..)) +
  geom_histogram(binwidth = 750, color = "white", fill = "grey75", size = 0.4, alpha = 1) +
  gg.opts + geom_density(col = "grey30", fill = "white", alpha = 0.4) +
  scale_x_continuous(label = comma, expand = c(0,0)) +
  scale_y_continuous(expand = c(0.03,0), breaks = seq(0, 0.15, 0.025)) +
  theme(axis.ticks = element_blank()) +
  bbc_style() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.ticks.x = ggplot2::element_line(color = "#cbcbcb"),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10))) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.2, 
               data = data.frame(x1 = eoo.alpha.ci[2], x2 = eoo.alpha.ci[3], y1 = 0.0, y2 = 0.0)) +
  geom_point(aes(x = x1, y = y1), colour = "black", size = 4,
             data = data.frame(x1 = eoo.alpha.ci[1], y1 = 0.0)) +
  ylab(expression(Density~("x"~10^3))) + 
  xlab("EOO")

fig.A2 <- eoo$all %>% 
  purrr:::map2_chr(.x = ., .y = aoo$all, .f = ~classify.threat(EOO = .x, AOO = .y)) %>% 
  table(.) %>% 
  as.data.frame(.) %>% 
  dplyr::rename(threat = ".", freq = "Freq") %>% 
  ggplot(data = ., aes(x = factor(threat), y = freq))+
  geom_bar(stat = "identity", width = 0.7, fill = "steelblue") +
  gg.opts + bbc_style() +
  xlab("Threat category") + ylab("Count") +
  theme(axis.text = element_text(size = 14, colour = "black")) +
  scale_y_continuous(breaks = seq(0, 900, 100), limits = c(0,900))

fig.A <- (fig.A1 | fig.A2) + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 18))

# # Save to file
# 
# ggsave(filename = "/Users/philippebouchet/Google Drive/Documents/projects/orcaella/AMMC/Papers/Distribution/FMS/Figures/Bouchet_etal_FigA.png", plot = fig.A, height = 4.5, width = 10)

# | Figure s2 ----------

pdf("out/Figure_S2a.pdf")
plot(kimb, col = "lightgrey", border = NA)
tracks[!names(tracks)%in%"OH_014"] %>% 
  purrr::map(.x = ., .f = ~ as(.x, "SpatialLines")) %>% 
  do.call(rbind, .) %>% 
  plot(., col = "#3FA0B3", add = TRUE)
dev.off()

pdf("out/Figure_S2b.pdf")
plot(kimb, col = "lightgrey", border = NA)
tracks["OH_014"] %>%
  purrr::map(.x = ., .f = ~ as(.x, "SpatialLines")) %>%
  do.call(rbind, .) %>%
  plot(., col = "black", add = TRUE)
dev.off()


# | Figure S3 ----------

# Coastline

bbo.sf <- kimb %>% sf::st_as_sf(.) %>% 
  sf::st_crop(x = ., xmin = 122.33, xmax = 122.38, ymin = -18, ymax = -17.97)

# Occurrence records

bbo.obs <- snub.onland %>% st_as_sf(x = ., coords = c("longitude", "latitude"), crs = CRSll) %>% 
  sf::st_crop(x = ., xmin = 122.327, xmax = 122.355, ymin = -18, ymax = -17.97)

# Buffers around points

bbo.buffer <- land.buffers %>% 
  sf::st_transform(x = ., crs = CRSll) %>% 
  sf::st_crop(x = ., xmin = 122.327, xmax = 122.355, ymin = -18, ymax = -17.97)

bbo.buffer <- bbo.buffer[2,]

# Plot

y_breaks <- seq(-17.995, -17.97, by = 0.005)
y_labels <- paste0(rev(seq(17.97, 17.995, by = 0.005)),"°S")
y_labels[c(2,4,6)] <- ""

bbo.map <- ggplot() + 
  geom_sf(data = bbo.buffer, fill = "lightblue", colour = "transparent") +
  geom_sf(data = bbo.sf) + 
  geom_sf(data = bbo.obs, colour = "deepskyblue3", size = 0.75) +
  geom_segment(aes(x = 122.3375, y = -17.97375, xend = 122.332, yend = -17.97375), size = 0.5, arrow = arrow(length = unit(0.05, "inches"), type = "closed")) +
  geom_segment(aes(x = 122.375, y = -17.9885, xend = 122.374, yend = -17.992), size = 0.5, arrow = arrow(length = unit(0.05, "inches"), type = "closed")) +
  geom_point(data = data.frame(longitude = 122.344233, latitude = -17.977067), aes(longitude, latitude), size = 2) +
  annotate("text", x = 122.335, y = -17.9725, label = "Broome") +
  annotate("text", x = 122.3443, y = -17.9755, label = "BBO", fontface = 'bold') +
  annotate("text", x = 122.37525, y = -17.9875, label = "Crab Creek") +
  annotate("text", x = 122.35, y = -17.9925, colour = "deepskyblue3", label = "Roebuck Bay", fontface = 'italic') +
  gg.opts + coord_sf(expand = FALSE, ) +
  xlab("") + ylab("") +
  theme(legend.key.height = unit(1.2, "cm"), 
        legend.title = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.border = element_rect(color = "black", size = 0.8, fill = "transparent")) +
  scale_y_continuous(breaks = y_breaks, 
                     labels = y_labels,
                     limits = c(-18, -17.97))

ggsave(filename = "out/FigS3.png", plot = bbo.map, device = "png", dpi = 650, width = 8, height = 5)

# | Figure s4 ----------

pdf("out/Figure_S4.pdf", height = 5, width = 9)
par(mfrow = c(1,2))
plot(density(snub.dat$depth), type = "l", main = "", xlab = "Depth (m)", ylab = "Density", col = "#45A6D6", lwd = 1.5)
text(-110, 0.08, expression(bold("A")), cex = 1.5)
plot(density(snub.dat$dfresh), type = "l", main = "", xlab = "Distance to river (km)", ylab = "Density", col = "#FFA908", lwd = 1.5)
text(130, 0.094, expression(bold("B")), cex = 1.5)
dev.off()