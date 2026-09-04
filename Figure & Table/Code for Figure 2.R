setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Mice/Male")
FDataA_M    <- read.csv(file = "FDataA_M.csv")

setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Mice/Female")
FDataA_F    <- read.csv(file = "FDataA_F.csv")

label_text_M <- paste("italic(R)^{2} == ", 0.88)
p1_M <- 
  ggplot(FDataA_M, aes(Log.OBS, Log.PRE)) + ## using log-sacle axis
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", linetype = "dashed",linewidth = 1, alpha = 0.8) +
  geom_point  (aes(shape   = as.factor(name)),size = 7, color = "#e3716e")  +
  annotation_logticks() +
  scale_y_continuous(limits = c(-8,8), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-8,8),labels = scales::math_format(10^.x))+
  coord_cartesian(clip = "off")
p1_M <- p1_M + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times New Roman"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text               = element_text (size   = 28, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 30, colour = "black", face = "bold"),   # label of axes
    axis.title.y = element_text(margin = margin(l = 10)),
    legend.position         =c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.title = element_text(size = 30, face = "bold"),  
    legend.text = element_text(size = 28)) +
  labs (title = "", x = "Observed values (ug/ml)",  y = "Predicted values (ug/ml or ug/g)", shape = "Sample")+
  annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = 1.2, vjust = 0.8, size = 15, family = "Times New Roman", colour = "black") +
  annotate("text", x = -4, y = 5, label = label_text_M , parse = TRUE, size = 15, color = "black", family = "Times New Roman")
print(p1_M)

label_text_F <- paste("italic(R)^{2} == ", 0.87)
p_F <- 
  ggplot(FDataA_F, aes(Log.OBS, Log.PRE)) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_point(aes(shape = as.factor(name), color = as.factor(name)), size = 7) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5),  labels = scales::math_format(10^.x)) +
  scale_color_manual(values = c("Kidney" = "steelblue4", "Lung" = "#7A70B5", "Liver" = "#FCBB44", "Plasma" = "#e3716e"), labels = c("Kidney", "Lung", "Liver", "Plasma")) +
  scale_shape_manual(values = c("Kidney" = 18, "Lung" = 17, "Liver" = 16, "Plasma" = 19), labels = c("Kidney", "Lung", "Liver", "Plasma")) +
  labs(title = "", x = "Observed values (ug/ml or ug/g)", y = "Predicted values (ug/ml or ug/g)", color = "Sample", shape = "Sample") +
  guides(color = guide_legend(title = "Sample"), shape = guide_legend(title = "Sample"))+
  coord_cartesian(clip = "off")
p1_F <- p_F + 
  theme(
    plot.background = element_rect(fill = "White"),
    text = element_text(family = "Times New Roman"),                                      # text front (Time new roman)
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "White"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 28, colour = "black", face = "bold"),       # tick labels along axes 
    axis.title = element_text(size = 30, colour = "black", face = "bold"),      # label of axes
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 28)
  ) +
  annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = 1.2, vjust = 0.8, size = 15, family = "Times New Roman", colour = "black") +
  annotate("text", x = -3, y = 4, label = label_text_F , parse = TRUE, size = 15, color = "black", family = "Times New Roman")
print(p1_F)

p1_F <- p1_F + theme(axis.title.y = element_blank())
p1_M <- p1_M + theme(plot.margin = margin(r = 20))
ggsave("Figure 2.tiff", grid.arrange(p1_M, p1_F, ncol = 2,  widths = c(1.02, 0.98)), path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 48, height = 24, units = "cm",dpi=320)
