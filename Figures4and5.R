library(ggplot2)

load(file = "/Users/delvaw/Documents/SimpactCyanExamples/figure4.ingredients.RData")
# revised.network.fortified, tree, dates, trans.and.nodes.long.enriched.combined

load(file = "/Users/delvaw/Documents/SimpactCyanExamples/figure5.ingredients.RData")
# input.output.F.long.revised.prev,
# input.output.F.long.revised.prev.model.average,
# input.output.F.long.revised.prev.average,
# darkcols,
# reds,
# blues,
# unaids.prev,
# input.output.F.long.revised.inc,
# input.output.F.long.revised.inc.model.average,
# input.output.F.long.revised.inc.average,
# shims2.df,
# input.output.F.long.revised.cov,
# input.output.F.long.revised.cov.model.average,
# input.output.F.long.revised.cov.average,
# unaids.art.cov

load(file = "/Users/delvaw/Documents/SimpactCyanExamples/fitGTRGI.top.RData")
# fitGTRGI.top, fitGTRGI.top.50.random, fitGTRGI.top.25.random, fitGTRGI.top.125.random
library(ape)
simseq.list <- fitGTRGI.top$data # This is the list of all simulated sequences
simseq.matrix <- as.character(simseq.list)
class(simseq.matrix) <- "matrix"
as.DNAbin(simseq.matrix)
dist.matrix <- dist.dna(as.DNAbin(as.character(simseq.list)), model = "raw", as.matrix = TRUE)

###
# FIGURE 4
###

cov.FaFc.revised.plot <- ggplot(data = input.output.F.long.revised.cov,
                                aes(year,
                                    art.cov,
                                    group = unique.id,
                                    colour = scenario.model)) +
  geom_line(data = input.output.F.long.revised.cov,
            linetype = "dashed",
            size = 0.3,
            show.legend = FALSE) +
  theme(legend.key = element_blank(),
        legend.background = element_blank()) +
  xlab("Time (years)") +
  ylab("ART coverage") +
  ylim(0, 1) +
  geom_line(data = input.output.F.long.revised.cov.model.average,
            aes(x = year,
                y = cov.av,
                group = scenario.model,
                col = scenario.model),
            linetype = "solid",
            size = 0.8) +
  geom_line(data = input.output.F.long.revised.cov.average,
            aes(x = year,
                y = cov.av,
                group = scenario,
                col = scenario.av),
            size = 1.5) +
  scale_color_manual(values = c("EAAA overall average" = darkcols[1],
                                "Counterfactual overall average" = darkcols[2],
                                "EAAA model 1" = reds[3],
                                "EAAA model 2" = reds[2],
                                "EAAA model 3" = reds[1],
                                "Counterfactual model 1" = blues[3],
                                "Counterfactual model 2" = blues[2],
                                "Counterfactual model 3" = blues[1]),
                     breaks = c("EAAA overall average",
                                "EAAA model 1",
                                "EAAA model 2",
                                "EAAA model 3",
                                "Counterfactual overall average",
                                "Counterfactual model 1",
                                "Counterfactual model 2",
                                "Counterfactual model 3"),
                     labels = c("EAAA overall average",
                                "EAAA model 1",
                                "EAAA model 2",
                                "EAAA model 3",
                                "Counterfactual overall average",
                                "Counterfactual model 1",
                                "Counterfactual model 2",
                                "Counterfactual model 3")) +
  geom_polygon(data = data.frame(x = c(unaids.art.cov$year+0.5,
                                       rev(unaids.art.cov$year+0.5)),
                                 y = c(unaids.art.cov$cov.lower/100,
                                       rev(unaids.art.cov$cov.upper/100))),
               inherit.aes = FALSE,
               aes(x = x, y = y),
               fill = "grey25",
               color = "grey25",
               alpha = 0.5,
               size = 1) +
  scale_x_continuous(breaks = seq(from = 2000,
                                  to = 2030,
                                  by = 5)) +
  theme(panel.background = element_rect(fill = "grey97"),
        axis.line.x = element_line(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

inc.FaFc.revised.plot <- ggplot(data = input.output.F.long.revised.inc,
                                aes(year, hivincidence, group = unique.id, colour = scenario.model)) +
  geom_line(linetype = "dashed",
            size = 0.3,
            show.legend = FALSE) +
  xlab("Time (years)") +
  ylab("HIV incidence") +
  geom_line(data = input.output.F.long.revised.inc.model.average,
            aes(x = year,
                y = inc.av,
                group = scenario.model,
                col = scenario.model),
            linetype = "solid",
            size = 0.8,
            show.legend = FALSE) +
  geom_line(data = input.output.F.long.revised.inc.average,
            aes(x = year,
                y = inc.av,
                group = scenario,
                col = scenario.av),
            size = 1.5) +
  scale_color_manual(values = c("EAAA overall average" = darkcols[1],
                                "Counterfactual overall average" = darkcols[2],
                                "EAAA model 1" = reds[3],
                                "EAAA model 2" = reds[2],
                                "EAAA model 3" = reds[1],
                                "Counterfactual model 1" = blues[3],
                                "Counterfactual model 2" = blues[2],
                                "Counterfactual model 3" = blues[1]),
                     breaks = c("EAAA overall average",
                                "EAAA model 1",
                                "EAAA model 2",
                                "EAAA model 3",
                                "Counterfactual overall average",
                                "Counterfactual model 1",
                                "Counterfactual model 2",
                                "Counterfactual model 3"),
                     labels = c("EAAA overall average",
                                "EAAA model 1",
                                "EAAA model 2",
                                "EAAA model 3",
                                "Counterfactual overall average",
                                "Counterfactual model 1",
                                "Counterfactual model 2",
                                "Counterfactual model 3")) +
  # overlaying the SHIMS II incidence estimate
  geom_polygon(data = shims2.df,
               inherit.aes = FALSE,
               aes(x = x, y = y),
               fill = "black",
               color = "black",
               alpha = 1,
               size = 4) +
  geom_segment(data = data.frame(x = 2014,
                                 xend = 2016.7,
                                 y = 0.0135,
                                 yend = 0.0135),
               inherit.aes = FALSE,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 1,
               arrow = arrow(type = "closed",
                             length = unit(0.03, "npc"))) +
  geom_label(data = shims2.df,
             inherit.aes = FALSE,
             aes(x = x, y = y, label = "UNAIDS estimate"),
             nudge_x = -11) +
  scale_x_continuous(breaks = seq(from = 1990,
                                  to = 2030,
                                  by = 5)) +
  theme(panel.background = element_rect(fill = "grey97"),
        axis.line.x = element_line(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

prev.FaFc.revised.plot <- ggplot(data = input.output.F.long.revised.prev,
                                 aes(year, hivprevalence, group = unique.id, colour = scenario.model)) +
  geom_line(linetype = "dashed",
            size = 0.3) +
  xlab("Time (years)") +
  ylab("HIV prevalence") +
  geom_line(data = input.output.F.long.revised.prev.model.average,
            aes(x = year,
                y = prev.av,
                group = scenario.model,
                col = scenario.model),
            linetype = "solid",
            size = 0.8,
            show.legend = FALSE) +
  geom_line(data = input.output.F.long.revised.prev.average,
            aes(x = year, y = prev.av, group = scenario, col = scenario.av),
            size = 1.5) +
  scale_color_manual(values = c("EAAA overall average" = darkcols[1],
                                "Counterfactual overall average" = darkcols[2],
                                "EAAA model 1" = reds[3],
                                "EAAA model 2" = reds[2],
                                "EAAA model 3" = reds[1],
                                "Counterfactual model 1" = blues[3],
                                "Counterfactual model 2" = blues[2],
                                "Counterfactual model 3" = blues[1]),
                     breaks = c("EAAA overall average",
                                "EAAA model 1",
                                "EAAA model 2",
                                "EAAA model 3",
                                "Counterfactual overall average",
                                "Counterfactual model 1",
                                "Counterfactual model 2",
                                "Counterfactual model 3"),
                     labels = c("EAAA overall average",
                                "EAAA model 1",
                                "EAAA model 2",
                                "EAAA model 3",
                                "Counterfactual overall average",
                                "Counterfactual model 1",
                                "Counterfactual model 2",
                                "Counterfactual model 3")) +
  geom_polygon(data = data.frame(x = c(unaids.prev$year+0.5,
                                       rev(unaids.prev$year+0.5)),
                                 y = c(unaids.prev$prev.lower/100,
                                       rev(unaids.prev$prev.upper/100))),
               inherit.aes = FALSE,
               aes(x = x, y = y),
               fill = "grey25",
               color = "grey25",
               alpha = 0.5,
               size = 1) +
  scale_x_continuous(breaks = seq(from = 1990,
                                  to = 2030,
                                  by = 5)) +
  theme(panel.background = element_rect(fill = "grey97"),
        axis.line.x = element_line(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")


###
# FIGURE 5
###
transmissionnetwork.plot <- ggplot(data = revised.network.fortified[2:nrow(revised.network.fortified), ]) +
  geomnet::geom_net(data = revised.network.fortified[2:nrow(revised.network.fortified), ],
                    aes(from_id = from_id,
                        to_id = to_id),
                    directed = TRUE,
                    size = 1, # 2.5
                    layout.alg = "kamadakawai",
                    alpha = 1,
                    arrowsize = 0.25, # 0.5
                    arrowgap = 0.004,
                    ecolour = "darkgrey",
                    colour = "black") +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "grey97")) +
  ylab("")


phylotree.plot <- ggtree::ggtree(tree,
                                 mrsd = dates[1],
                                 size = 0.05) + 
  ggtree::theme_tree2() +
  theme_grey() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "grey97")) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) +
  xlab("Time") +
  ylab("")


transandnodescombined.plot <- ggplot(data = dplyr::filter(trans.and.nodes.long.enriched.combined,
                                                          calendaryear >= 1985,
                                                          calendaryear <= 2017),
                                     aes(x = calendaryear,
                                         y = percentage.events,
                                         colour = Scenario)) +
  geom_point() +
  geom_line(size = 1,
            alpha = 0.75) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time") +
  ylab("Probability density") +
  theme(axis.line.x = element_line(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "grey97")) +
  scale_color_manual(values = c("Transmission events" = custom.colours[1],
                                "100% sequence coverage" = custom.colours[2],
                                "Random 50% sequence coverage" = custom.colours[3],
                                "Random 25% sequence coverage" = custom.colours[4],
                                "Random 12.5% sequence coverage" = custom.colours[5])) 


###
# Combined panels
###
Figure4 <- ggpubr::ggarrange(cov.FaFc.revised.plot,
                  inc.FaFc.revised.plot,
                  prev.FaFc.revised.plot, 
          labels = c("a", "b", "c"),
          nrow = 1, ncol = 3)
ggplot2::ggsave(filename = "Figure4.pdf",
       plot = Figure4,
       path = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision2/plots",
       width = 36, height = 10, units = "cm")

Figure5 <- ggpubr::ggarrange(transmissionnetwork.plot,
                             phylotree.plot,
                             transandnodescombined.plot,
                             labels = c("a", "b", "c"),
                             nrow = 1, ncol = 3,
                             widths = c(1, 1, 1.5))
ggplot2::ggsave(filename = "Figure5.pdf",
                plot = Figure5,
                path = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/plots",
                width = 36, height = 10, units = "cm")

# Calibration output Table (for second round of revisions)
load(file = "/Users/delvaw/Documents/SimpactCyanExamples/figure5.ingredients.RData")

load(file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/SimpactPaperEAAAExample.Revision2.RData")
# This loads the "input.output.F.wide.revised" dataframe

# The vector of target features was:
features.pop.growth <- exp(0.015)
features.hiv.prev <- c(0.143, 0.008, 0.315, 0.066, 0.467, 0.213, 0.538, 0.366, 0.491, 0.47, 0.397, 0.455, 0.316, 0.425)
features.hiv.inc <- exp(c(0.038, 0.008, 0.043, 0.016, 0.02, 0.026, 0.027, 0.031, 0.04, 0.004, 0.021, 0.012, 0.012, 0))
features.art.cov <- c(0.37, 0.40, 0.44, 0.49, 0.58, 0.67, 0.76, 0.85)
features.vl.suppr <- 0.74
features.unaids.prev <- c(0.017, 0.033, 0.057, 0.088, 0.124, 0.161, 0.194, 0.220, 0.239, 0.251, 0.258, 0.261, 0.261, 0.259, 0.257, 0.255, 0.256, 0.259, 0.263, 0.268, 0.274, 0.278, 0.282, 0.284, 0.284, 0.283, 0.279, 0.274)
target.features.EAAA.for.table <- c(log(features.pop.growth), features.hiv.prev, round(log(features.hiv.inc), 3), features.art.cov, features.vl.suppr, features.unaids.prev)

# # MODEL OUTPUT:
# 1 pop growth +
# 14 age-gender-specific prev +
# 14 age-gender-specific inc +
# 65 ART coverage +       # 20:52, by 0.5
# 1 vl suppression +
# 180 incidence estimates + # from = 7.25, to = 52, by = 0.25
# 180 incident cases +
# 129 art cases +
# 1 SHIMS1 prev +
# 1 SHIMS1 inc +
# 1 SHIMS2 inc +
# 92 annual HIV prev     # 6.5:52, by = 0.5
# TOTAL: 679 output statistics

input.output.F.wide.revised.for.prev <- input.output.F.wide.revised[, c((x.offset+30+65+1+360+129+3):(x.offset+(30+65+360+129+3+92)), (index.end - 2):index.end)] # Making a copy for prevalence plot


only.output.F.wide.revised <- input.output.F.wide.revised[, c(21:699, 701:703)]

# The columns we need are:
# 1 (pop growth)
# 2:29 (age-gender specific prev and inc)
# 51:2:65  # 30~1 (22:2:36 of the 65 measurements)
# 95 (vl suppression)
# 596:2:650 # We need 1990.5 until 2017.5
output.for.calib.table <- only.output.F.wide.revised[only.output.F.wide.revised$scenario == "EAAA ", c(1,
                                                                                                       2:29,
                                                                                                       seq(51,65,by=2),
                                                                                                       95,
                                                                                                       seq(596,c(596+2*27),by=2),
                                                                                                       682)]
output.for.calib.table[, 1] <- log(output.for.calib.table[, 1])
output.for.calib.table[, 16:29] <- log(output.for.calib.table[, 16:29])

library(tidyverse)
by_model <- output.for.calib.table %>%
  group_by(model)

calib.tibble <- by_model %>%
  summarise_all(mean) %>%
  t() %>%
  cbind(., c("target feature", target.features.EAAA.for.table)) 

calib.df <- as.data.frame(calib.tibble[2:67, ])
indx <- sapply(calib.df, is.factor)
calib.df[indx] <- lapply(calib.df[indx], function(x) as.numeric(as.character(x)))
calib.df <-  round(calib.df, 3) / 100

write.csv2(calib.df, file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision2/calib.csv")



