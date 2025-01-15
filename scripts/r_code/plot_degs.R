#!/usr/bin/env Rscript

library(ggtext)

# prepare data for visualization
p_data <- as_tibble(diff_exp_genes) %>%
  group_by(timepoint, regulation) %>%
  reframe(num_genes = n()) %>%
  mutate(pos_shift = ifelse(regulation == "up", 0.2, -0.2),
         timepoint = factor(timepoint, levels = c("t4", "t12", "t24", "t72")))

# plot of differentially expressed gene tallies
deg_bar_plt <- p_data %>% 
  ggplot(aes(timepoint, num_genes, fill = regulation)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = num_genes), position = position_dodge(width = 0.9),
            vjust = -0.5, fontface = "bold", size = 8) +
  scale_fill_manual(values = c("#0000ff", "#ff0000"),
                    breaks = c("down", "up"),
                    labels = c("Downregulated", "Upregulated")) +
  scale_y_continuous(limits = c(0, 2000),
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.16, 0.16),
                   breaks = c("t4", "t12", "t24", "t72"),
                   labels = c("4-hpi", "12-hpi", "24-hpi", "72-hpi")) +
  labs(title = "Differentially Expressed Genes of<br> THEV-infected B-cells",
       y = "Number of Genes",
       x = element_blank(),
       fill = element_blank()) +
  theme_classic() +
  theme(plot.margin = margin(rep(20, 4)),
        plot.title = element_markdown(face = 'bold',
                                      size = 23,
                                      colour = 'black',
                                      hjust = 0,
                                      lineheight = -10,
                                      vjust = 0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major.y = element_line(color = 'grey50',
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.minor = element_line(linewidth = 0.1, colour = "grey50",
                                        linetype = "dashed"),
        #adjust axis
        axis.text.x = element_markdown(size = 18, colour = 'black', face = 'bold',
                                       margin = margin(t = 10)),
        axis.text.y = element_text(size = 18, colour = 'black', face = 'bold',
                                   margin = margin(l = 10)),
        axis.title.y = element_text(size = 22,
                                    face = 'bold',
                                    color = 'black'),
        axis.ticks = element_blank(),
        legend.justification = c(0,1),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.95),
        legend.box.background = element_rect(color = "grey40"),
        legend.text = element_text(face = 'bold', 
                                   size = 18,
                                   family = 'Arial'),
        legend.background = element_blank(),
        legend.margin = margin(rep(10, 4)),
        legend.key.size = unit(1.2, "cm")
  )
