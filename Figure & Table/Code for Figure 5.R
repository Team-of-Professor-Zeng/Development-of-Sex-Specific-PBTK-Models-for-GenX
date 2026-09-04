library(readxl)
library(ggplot2)
library(ggtext)
library(dplyr)
library(gridExtra)
library(sysfonts)
library(showtext)
font_add("Times New Roman", regular = "C:/Windows/Fonts/times.ttf",
         bold = "C:/Windows/Fonts/timesbd.ttf",
         italic = "C:/Windows/Fonts/timesi.ttf",
         bolditalic = "C:/Windows/Fonts/timesbi.ttf")
showtext_auto()

setwd("C:/Users/15960/Desktop/2025GenX/Modfit/Human")
HQ  <- read.csv(file = "HQ-Plot.csv")
#####################################################################################################################################
## ASC
HQ1 <- HQ %>% filter(Type == "ASC" & Sex == "Male"   & Exposure == "Chronic"   )  %>% select(Regions = "Regions", Median = "Median")
HQ3 <- HQ %>% filter(Type == "ASC" & Sex == "Female" & Exposure == "Chronic"   )  %>% select(Regions = "Regions", Median = "Median")
HQ1$Regions <- factor(HQ1$Regions,
levels = c("Dong et al. (2023), Shandong (Huantai)", "Feng et al. (2021), Shandong (Huantai)",
"Dong et al. (2024), China (Nationwide)", "Dong et al. (2024), China (North Coast)",
"Dong et al. (2024), China (Northeast)", "Dong et al. (2024), China (Northwest)",
"Dong et al. (2024), China (Midland)", "Dong et al. (2024), China (South Coast)",
"Xu et al. (2021), Zhejiang (Hangzhou)"))
# Define regions requiring a superscript
special_regions <- c("Dong et al. (2023), Shandong (Huantai)", "Feng et al. (2021), Shandong (Huantai)")
# Male uses log scale; set a positive X-axis lower limit to avoid negative values
min_val1 <- min(HQ1$Median) * 0.9
max_val1 <- max(HQ1$Median) * 1.15
x_limits1 <- c(min_val1, max_val1)
p1 <- ggplot(HQ1, aes(Median, Regions)) +
geom_vline(xintercept = 1, linetype = 2, color = "black", linewidth = 1) +
annotate("rect", xmin = x_limits1[1], xmax = 1, ymin = -Inf, ymax = Inf, fill = "#B8DBB3", alpha = 0.2) +
annotate("rect", xmin = 1, xmax = x_limits1[2], ymin = -Inf, ymax = Inf, fill = "#F7A6AC", alpha = 0.2) +
scale_x_log10(limits = x_limits1, expand = expansion(mult = c(0, 0)),
              breaks = c(0.01, 0.1, 1, max(HQ1$Median)),
              labels = c("0.01", "0.1", "1", sprintf("%.2f", max(HQ1$Median)))) +
scale_y_discrete(limits = rev,
labels = function(x) ifelse(x %in% special_regions, paste0(x, "<sup>a</sup>"), x),
expand = c(0.12, 0)) +
geom_point(aes(color = Median > 1, shape = Median > 1), size = 8) +
scale_color_manual(values = c("TRUE" = "#e3716e", "FALSE" = "#0A8B4F")) +
scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 18)) +
annotate("text", x = 0.2, y = Inf, label = "Negligible Risk", hjust = 0.5, vjust = -0.6, size = 9,
colour = "black", family = "Times New Roman", fontface = "bold") +
annotate("text", x = 2.6, y = Inf, label = "Potential Risk", hjust = 0.5, vjust = -0.6, size = 9,
colour = "black", family = "Times New Roman", fontface = "bold") +
annotate("text", x = x_limits1[1], y = Inf, label = "(A) Male",
hjust = 3.5, vjust = 0.2, size = 14, colour = "black", family = "Times New Roman", fontface = "bold") +
coord_cartesian(clip = "off")
p1 <- p1 +
theme(
plot.background  = element_rect(fill = "transparent"),
text             = element_text(family = "Times New Roman"),
panel.border     = element_rect(colour = "black", fill = NA, linewidth = 2),
panel.background = element_rect(fill = "White"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text        = element_text(size = 25, colour = "black", face = "bold"),
axis.text.y      = element_markdown(size = 27, margin = margin(r = 40), face = "plain", hjust = 0),
axis.title.y     = element_text(margin = margin(r = 30)),
axis.title       = element_text(size = 30, colour = "black", face = "bold"),
plot.margin      = margin(t = 40, r = 40, b = 10, l = 80, unit = "pt"), 
legend.position  = "none"
) +
labs(x = "Hazard Quotient (HQ)", y = "Regions")
print(p1)
#####################################################################
## ASC - Female
HQ3 <- HQ %>% filter(Type == "ASC" & Sex == "Female" & Exposure == "Chronic") %>% select(Regions, Median)
HQ3$Regions <- factor(HQ3$Regions,
levels = c("Dong et al. (2023), Shandong (Huantai)", "Feng et al. (2021), Shandong (Huantai)",
"Dong et al. (2024), China (Nationwide)", "Dong et al. (2024), China (North Coast)",
"Dong et al. (2024), China (Northeast)", "Dong et al. (2024), China (Northwest)",
"Dong et al. (2024), China (Midland)", "Dong et al. (2024), China (South Coast)",
"Xu et al. (2021), Zhejiang (Hangzhou)"))
max_val3 <- max(abs(HQ3$Median - 1)) * 1.05
x_limits3 <- c(1 - max_val3, 1 + max_val3)
p3 <- ggplot(HQ3, aes(Median, Regions)) +
geom_vline(xintercept = 1, linetype = 2, color = "black", linewidth = 1) +
annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#B8DBB3", alpha = 0.2) +
annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F7A6AC", alpha = 0.2) +
scale_x_continuous(limits = x_limits3, expand = c(0.1, 0)) +
scale_y_discrete(limits = rev,
labels = function(x) ifelse(x %in% special_regions, paste0(x, "<sup>a</sup>"), x),
expand = c(0.12, 0)) +
geom_point(aes(color = Median > 1, shape = Median > 1), size = 8) +
scale_color_manual(values = c("TRUE" = "#e3716e", "FALSE" = "#0A8B4F")) +
scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 18)) +
annotate("text", x = 0.4, y = Inf, label = "Negligible Risk", hjust = 0.5, vjust = -0.6, size = 9,
colour = "black", family = "Times New Roman", fontface = "bold") +
annotate("text", x = 1.6, y = Inf, label = "Potential Risk", hjust = 0.5, vjust = -0.6, size = 9,
colour = "black", family = "Times New Roman", fontface = "bold") +
annotate("text", x = -Inf, y = Inf, label = "(B) Female",
hjust = 2.8, vjust = 0.2, size = 14, colour = "black", family = "Times New Roman", fontface = "bold") +
coord_cartesian(clip = "off")
p3 <- p3 +
theme(
plot.background  = element_rect(fill = "transparent"),
text             = element_text(family = "Times New Roman"),
panel.border     = element_rect(colour = "black", fill = NA, linewidth = 2),
panel.background = element_rect(fill = "White"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text        = element_text(size = 25, colour = "black", face = "bold"),
axis.text.y      = element_markdown(size = 27, margin = margin(r = 40), face = "plain", hjust = 0),
axis.title.y     = element_text(margin = margin(r = 30)),
axis.title       = element_text(size = 30, colour = "black", face = "bold"),
plot.margin      = margin(t = 40, r = 40, b = 10, l = 80, unit = "pt"),
legend.position  = "none"
) +
labs(x = "Hazard Quotient (HQ)", y = "Regions")
print(p3)
#####################################################################
## Save figure (ensure directory exists)
art_path <- "C:/Users/15960/Desktop/2025GenX/Artwork"
if(!dir.exists(art_path)) dir.create(art_path, recursive = TRUE)
ggsave("Figure_5_grey.tiff",
plot = grid.arrange(p1, p3, nrow = 2),
path = art_path,
width = 40, height = 36, units = "cm", dpi = 320)
