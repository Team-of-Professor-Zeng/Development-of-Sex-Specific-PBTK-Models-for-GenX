DataA_M_combined <- FDataA_M %>% mutate(Sex = "Male")
FDataA_F_combined <- FDataA_F %>% mutate(Sex = "Female")
FDataA_combined <- bind_rows(FDataA_M_combined, FDataA_F_combined)

pct_data <- FDataA_combined %>%
  group_by(Sex) %>%
  summarize(
    pct_2e = mean(OPR >= 1/2 & OPR <= 2) * 100,
    pct_3e = mean(OPR >= 1/3 & OPR <= 3) * 100
  )

p_combined <- ggplot(FDataA_combined, aes(x = OPR, fill = Sex)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), 
                 bins = 30, color = "black", alpha = 0.7, 
                 position = "identity") +  
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0,0)) +
  scale_fill_manual(values = c("Male" = "steelblue", "Female" = "#e3716e")) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e3716e", linewidth = 1) + 
  geom_vline(xintercept = c(1/2, 2), linetype = "dashed", color = "black", linewidth = 0.8) + 
  geom_vline(xintercept = c(1/3, 3), linetype = "dashed", color = "grey", linewidth = 0.8) + 
  geom_text(data = pct_data, 
            aes(x = 0.005, y = 32, 
                label = paste0("% 2e: ", round(pct_2e, 1), "%")),
            hjust = 0, vjust = 1, size = 6, color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  geom_text(data = pct_data, 
            aes(x = 0.005, y = 30, 
                label = paste0("% 3e: ", round(pct_3e, 1), "%")),
            hjust = 0, vjust = 1, size = 6, color = "grey", fontface = "bold",
            inherit.aes = FALSE) +
  labs(x = "Predicted / Observed", y = "Proportion (%)", fill = "Sex") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.background = element_rect(fill = "white"),
    text = element_text(family = "Times New Roman"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 20, colour = "black", face = "bold"),
    axis.title = element_text(size = 22, colour = "black", face = "bold"),
    legend.position = c(0.95, 0.95),  
    legend.justification = c(1, 1),   
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    legend.text = element_text(size = 14, face = "bold",
                               margin = margin(t = 8, b = 8, unit = "pt")),  
    legend.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.8, 0.6, 0.3, 0.5, "cm"))
print(p_combined)


pct_data <- FDataA_combined %>%
  group_by(Sex) %>%
  summarize(
    pct_2e = mean(OPR >= 1/2 & OPR <= 2) * 100,
    pct_3e = mean(OPR >= 1/3 & OPR <= 3) * 100
  )

legend_data <- data.frame(
  Sex = c("Male", "Female"),
  pct_2e = pct_data$pct_2e,
  pct_3e = pct_data$pct_3e
)

p_combined <- ggplot(FDataA_combined, aes(x = OPR, fill = Sex)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), 
                 bins = 30, color = "black", alpha = 0.7, 
                 position = "identity") +  
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0,0)) +
  scale_fill_manual(values = c("Male" = "steelblue", "Female" = "#e3716e")) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e3716e", linewidth = 1) + 
  geom_vline(xintercept = c(1/2, 2), linetype = "dashed", color = "black", linewidth = 0.8) + 
  geom_vline(xintercept = c(1/3, 3), linetype = "dashed", color = "grey", linewidth = 0.8) + 
  labs(x = "Predicted / Observed", y = "Proportion (%)") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.background = element_rect(fill = "white"),
    text = element_text(family = "Times New Roman"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 20, colour = "black", face = "bold"),
    axis.title = element_text(size = 22, colour = "black", face = "bold"),
    legend.position = "none", 
    plot.margin = margin(0.8, 0.6, 0.3, 0.5, "cm"))

custom_legend <- ggplot(legend_data, aes(x = 1, y = Sex)) +
  geom_tile(aes(fill = Sex), width = 0.1, height = 0.3, color = "black") +
  geom_text(aes(label = Sex),
            hjust = 0, vjust = 0.5, size = 6, fontface = "bold", color = "black",
            position = position_nudge(x = 0.1)) +
  geom_text(aes(label = paste0("\n%2e: ", round(pct_2e, 1), "%\n%3e: ", round(pct_3e, 1), "%")),
            hjust = -0.05, vjust = 0.9, size = 5, fontface = "bold", color = "black",
            position = position_nudge(x = 0.1)) +
  annotate("text", x = 0.97, y = 2.47, label = "Sex", 
           size = 7, fontface = "bold", hjust = 0.12) +
  scale_fill_manual(values = c("Male" = "steelblue", "Female" = "#e3716e")) +
  scale_y_discrete(limits = rev) +  
  coord_cartesian(xlim = c(0.8, 2.5), ylim = c(0.5, 2.5)) + 
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "NA", color = "NA", linewidth = 1),
    plot.margin = margin(2, 5, 2, 5)
  ) 

library(patchwork)
final_plot <- p_combined + 
  inset_element(custom_legend, 
                left = -0.05, bottom = 0.6, right = 0.8, top = 1.05) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white")))

print(final_plot)

ggsave("Figure S1.tiff", plot = final_plot,
       path = "C:/Users/15960/Desktop/2025GenX/Artwork",
       width = 20,   height = 20,  units = "cm", dpi = 320)