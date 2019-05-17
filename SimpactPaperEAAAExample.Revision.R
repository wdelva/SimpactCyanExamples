setwd("/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/")
# path/to/your/working_directory/") # Change this to your own working directory
library(devtools)
library(RSimpactCyan)
install_github("wdelva/RSimpactHelp")
library(RSimpactHelper)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(metafolio)

# Loading the output from the Lenormand calibration
EAAA.SciRep.revision.Seq.loaded.object <- load(file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/EAAA.SciRep.revision2.Seq.RData")

Seq.object <- EAAA.SciRep.revision2.Seq

# The vector of target features was:
features.pop.growth <- exp(0.015)
features.hiv.prev <- c(0.143, 0.008, 0.315, 0.066, 0.467, 0.213, 0.538, 0.366, 0.491, 0.47, 0.397, 0.455, 0.316, 0.425)
features.hiv.inc <- exp(c(0.038, 0.008, 0.043, 0.016, 0.02, 0.026, 0.027, 0.031, 0.04, 0.004, 0.021, 0.012, 0.012, 0))
features.art.cov <- c(0.37, 0.40, 0.44, 0.49, 0.58, 0.67, 0.76, 0.85)
features.vl.suppr <- 0.74
features.unaids.prev <- c(0.017, 0.033, 0.057, 0.088, 0.124, 0.161, 0.194, 0.220, 0.239, 0.251, 0.258, 0.261, 0.261, 0.259, 0.257, 0.255, 0.256, 0.259, 0.263, 0.268, 0.274, 0.278, 0.282, 0.284, 0.284, 0.283, 0.279, 0.274)
target.features.EAAA <- c(features.pop.growth, features.hiv.prev, features.hiv.inc, features.art.cov, features.vl.suppr, features.unaids.prev)

# To find the best fitting model, we compute the RMSE for the 250 parameter combinations of the posterior.
posterior.size <- nrow(Seq.object$stats)
RMSE.vect <- rep(NA, posterior.size)
for (i.posterior in 1:posterior.size){
  fitting.features.EAAA <- Seq.object$stats[i.posterior, ]
  # The RMSE compared to the target statistics was:
  RMSE.vect[i.posterior] <- sqrt(sum(((fitting.features.EAAA - target.features.EAAA)/target.features.EAAA)^2) / length(target.features.EAAA))
}

index.bestfit <- which(RMSE.vect == min(RMSE.vect))
bestfitting.features.EAAA <- Seq.object$stats[index.bestfit, ]

# Visualising the relative distance:
# hist((bestfitting.features.EAAA - target.features.EAAA)/target.features.EAAA, 50)
# plot(1:66, (bestfitting.features.EAAA - target.features.EAAA)/target.features.EAAA)

# Using the best-fitting model from the calibration
inputvector.EAAA.example <- Seq.object$param[index.bestfit, ]
# Or we could use multiple parameter combinations, e.g. the 3 smallest RMSEs
indices.best3fits <- order(RMSE.vect)[1:3]
RMSE.vect[indices.best3fits]
# 0.08598585 0.10758438 0.13292604  These are the RMSrelE for the best 3 fitting models
inputvectors.EAAA.example <- Seq.object$param[indices.best3fits, ] #[indices.best10fits, ]

# The original submission used this input vector:
# inputvector.EAAA.example <- c(1.1,
#                               0.4,
#                               1,
#                               3.5,
#                               0.8,
#                               0.3,
#                               45,
#                               30,
#                               -1,
#                               3.8,
#                               -0.5,
#                               -0.45,
#                               -2.5,
#                               -1,
#                               -5,
#                               2,
#                               2,
#                               1,
#                               0)

reps <- 10    # Original submission used 20 repititions

# The design matrix contains multiple models (i.e. parameter combinations)
EAAA.multimodel.design.matrix <- matrix(rep(inputvectors.EAAA.example,
                                 each = reps),
                             nrow = reps * nrow(inputvectors.EAAA.example),
                             byrow = FALSE)

EAAA.revised.F.a.output <- simpact.parallel(model = EAAA.revised.F.a.wrapper,
                                            actual.input.matrix = EAAA.multimodel.design.matrix,
                                            seed_count = 0,
                                            n_cluster = 4)

EAAA.revised.F.c.output <- simpact.parallel(model = EAAA.revised.F.c.wrapper,
                                            actual.input.matrix = EAAA.multimodel.design.matrix,
                                            seed_count = 0,
                                            n_cluster = 4)

input.df <- data.frame(seed = 1:nrow(EAAA.multimodel.design.matrix),
                       EAAA.multimodel.design.matrix)

input.output.Fa.df <- cbind(input.df, as.data.frame(EAAA.revised.F.a.output))
input.output.Fc.df <- cbind(input.df, as.data.frame(EAAA.revised.F.c.output))
input.output.Fa.df$scenario <- "EAAA "
input.output.Fc.df$scenario <- "Counterfactual "


## combining Fa and Fc datasets
input.output.F.wide.revised <- rbind(input.output.Fa.df,
                             input.output.Fc.df)
input.output.F.wide.revised$unique.id <- 1:nrow(input.output.F.wide.revised)

# Grouping variable to group simulations that used the same input parameters
input.output.F.wide.revised$model <- rep(paste0("model ", 1:nrow(inputvectors.EAAA.example)),
   each = reps)

x.offset <- 20 # seed + 19

# MODEL OUTPUT:
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

# EAAA.revised.F.a.output and EAAA.revised.F.c.output have 680 columns because simpact.parallel adds the seedcount as well.
# Before conducting the actual analysis, we need to filter out any runs for which in 1 or both of the scenarios, NA values exist.

Fa.complete.rows <- complete.cases(input.output.Fa.df)
Fc.complete.rows <- complete.cases(input.output.Fc.df)

complete.runs <- Fa.complete.rows & Fc.complete.rows

input.output.F.wide.revised <- input.output.F.wide.revised[c(complete.runs, complete.runs), ]

save(input.output.F.wide.revised, file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/SimpactPaperEAAAExample.Revision2.RData")
revision2.object <- load(file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/SimpactPaperEAAAExample.Revision2.RData")


# Colour scheme
# dark red #7D0E0F

# light red #C90001
# pinkish red #E85C5D
# baby pink #FFB7F8


# dark blue #00427A # #20496B

# light blue #096CBF
# lighter blue #4FAFFF
# baby blue #A3D5FF

cols <- c("#E41A1C", "#377EB8")
darkcols <- c("#7D0E0F", "#00427A")
reds <- c("#C90001", "#E85C5D", "#FFB7F8")
blues <- c("#096CBF", "#4FAFFF", "#A3D5FF")

index.end <- ncol(input.output.F.wide.revised)

#### ART COVERAGE PLOT
input.output.F.wide.revised.for.cov <- input.output.F.wide.revised[, c((x.offset+30):(x.offset+(30+64)), (index.end - 2):index.end)] # Making a copy for ART coverage plot
names(input.output.F.wide.revised.for.cov)[1:65] <- seq(from = 2000,
                                                to = 2032,
                                                by = 0.5)

input.output.F.long.revised.cov <- gather(data = input.output.F.wide.revised.for.cov,
                                  key = year,
                                  value = art.cov,
                                  1:65)
input.output.F.long.revised.cov$year <- as.numeric(input.output.F.long.revised.cov$year)
input.output.F.long.revised.cov$scenario.model <- paste0(input.output.F.long.revised.cov$scenario,
                                                         input.output.F.long.revised.cov$model)

input.output.F.long.revised.cov$scenario.model <- factor(input.output.F.long.revised.cov$scenario.model,
                                                         levels = c("EAAA model 1",
                                                                    "EAAA model 2",
                                                                    "EAAA model 3",
                                                                    "Counterfactual model 1",
                                                                    "Counterfactual model 2",
                                                                    "Counterfactual model 3"))

## Adding average ART coverage per year and scenario, averaged over reps
input.output.F.long.revised.cov.average <- input.output.F.long.revised.cov %>%
  group_by(., scenario, year) %>%
  summarise(cov.av = mean(art.cov))
input.output.F.long.revised.cov.average$scenario.av <- paste0(input.output.F.long.revised.cov.average$scenario,
                                                              "overall average")

## Adding average ART coverage per year, scenario, and model, averaged over reps
input.output.F.long.revised.cov.model.average <- input.output.F.long.revised.cov %>%
  group_by(., scenario, model, year) %>%
  summarise(cov.av = mean(art.cov))
input.output.F.long.revised.cov.model.average$scenario.model <- paste0(input.output.F.long.revised.cov.model.average$scenario,
                                                                 input.output.F.long.revised.cov.model.average$model)


input.output.F.long.revised.cov.model.average$scenario.model <- factor(input.output.F.long.revised.cov.model.average$scenario.model,
                                                          levels = c("EAAA model 1",
                                                                         "EAAA model 2",
                                                                         "EAAA model 3",
                                                                         "Counterfactual model 1",
                                                                         "Counterfactual model 2",
                                                                         "Counterfactual model 3"))

input.output.F.long.revised.cov.average$scenario.av <- factor(input.output.F.long.revised.cov.average$scenario.av,
                                                              levels = c("EAAA overall average",
                                                                         "Counterfactual overall average"))

# FaFc ART coverage plot

# Overlaying the UNAIDS ART coverage estimates
setwd("/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/")
unaids.art.cov <- read.csv(file = "Eswatini_UNAIDS_2018_Treatment cascade_National.csv",
                           sep = ";",
                           skip = 2,
                           header = FALSE,
                           col.names = c("year",
                                         "cov.est",
                                         "cov.lower",
                                         "cov.upper"))
unaids.art.cov$unique.id <- NA
unaids.art.cov$scenario <- NA


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
  theme(panel.background = element_rect(fill = "grey97"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
plot(cov.FaFc.revised.plot)

# Saving the ART coverage plot
ggsave(filename = "cov_FaFc_3models_revised_plot.pdf",
       plot = cov.FaFc.revised.plot,
       path = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/plots",
       width = 16, height = 10, units = "cm")


#### HIV INCIDENCE PLOT
input.output.F.wide.revised.for.inc <- input.output.F.wide.revised[, c((x.offset+30+65+1):(x.offset+(30+65+180)), (index.end - 2):index.end)] # Making a copy for incidence plot
# from = 7.25, to = 52, by = 0.25
names(input.output.F.wide.revised.for.inc)[1:180] <- seq(from = 1987.25, to = 2032, by = 0.25)

input.output.F.long.revised.inc <- gather(data = input.output.F.wide.revised.for.inc,
                                  key = year,
                                  value = hivincidence,
                                  1:180)
input.output.F.long.revised.inc$year <- as.numeric(input.output.F.long.revised.inc$year)
input.output.F.long.revised.inc$scenario.model <- paste0(input.output.F.long.revised.inc$scenario,
                                                         input.output.F.long.revised.inc$model)

## Adding average HIV incidence per time point and scenario, averaged over reps
input.output.F.long.revised.inc.average <- input.output.F.long.revised.inc %>%
  group_by(., scenario, year) %>%
  summarise(inc.av = mean(hivincidence))
input.output.F.long.revised.inc.average$scenario.av <- paste0(input.output.F.long.revised.inc.average$scenario, "overall average")

## Adding average HIV incidence per time point, scenario, and model, averaged over reps
input.output.F.long.revised.inc.model.average <- input.output.F.long.revised.inc %>%
  group_by(., scenario, model, year) %>%
  summarise(inc.av = mean(hivincidence))
input.output.F.long.revised.inc.model.average$scenario.model <- paste0(input.output.F.long.revised.inc.model.average$scenario,
                                                                       input.output.F.long.revised.inc.model.average$model)


# FaFc incidence plot
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
  # overlaying the SHIMS II incidence estimate (1.48 [95% CI: 0.96 - 2.00] )
  geom_polygon(data = data.frame(x = c(2017.5, 2017.5, 2017.5, 2017.5),
                                 y = c(1.35, 1.35, 1.35, 1.35)/100),
               inherit.aes = FALSE,
               aes(x = x, y = y),
               fill = "black",
               color = "black",
               alpha = 1,
               size = 4) +
  theme(panel.background = element_rect(fill = "grey97"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
plot(inc.FaFc.revised.plot)
# Saving the incidence plot
ggsave(filename = "inc_FaFc_3models_revised_plot.pdf",
       plot = inc.FaFc.revised.plot,
       path = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/plots",
       width = 16, height = 10, units = "cm")


# SUMMARISING IMPACT OF EAAA ON HIV INCIDENCE 
# If we turn the long into a wide dataset, with hivincidence for Fa and hivincidence for Fc,
# we can easily subtract the one from the other.

input.output.F.long.revised.inc$run.id <- input.output.F.long.revised.inc$unique.id %% reps

scenariowide.model.inc <- dplyr::select(input.output.F.long.revised.inc,
                                                               -c(scenario,
                                                                  unique.id,
                                                                  model))




impact.summary.df <- tidyr::spread(scenariowide.model.inc,
                                   key = scenario.model,
                                   value = hivincidence)
impact.summary.df <- dplyr::rename(impact.summary.df,
                                   Fa.model1 = "EAAA model 1",
                                   Fa.model2 = "EAAA model 2",
                                   Fa.model3 = "EAAA model 3",
                                   Fc.model1 = "Counterfactual model 1",
                                   Fc.model2 = "Counterfactual model 2",
                                   Fc.model3 = "Counterfactual model 3")


impact.summary.df <- dplyr::mutate(impact.summary.df,
         relreduc1 = 1 - Fa.model1/Fc.model1,
         relreduc2 = 1 - Fa.model2/Fc.model2,
         relreduc3 = 1 - Fa.model3/Fc.model3)

impact.summary.df.long <- tidyr::gather(impact.summary.df,
                                        key = model,
                                        value = relreduc,
                                        relreduc1:relreduc3)
time.summary.grouped <- group_by(impact.summary.df.long,
                         year, model) %>%
  summarise(mean.relreduc = mean(relreduc))

summary(time.summary.grouped$mean.relreduc[time.summary.grouped$year>=2032.0])


impact.plot <- ggplot(data = dplyr::filter(time.summary.grouped,
                                           year >= 2016.75),
                       aes(x = year, y = mean.relreduc, group = model)) +
  geom_line(aes(colour = model))
plot(impact.plot)


## Model-averaged impact of EAAA

time.modelav.summary.grouped <- group_by(time.summary.grouped,
                                 year) %>%
  summarise(model.averaged.relreduc = mean(mean.relreduc))

impact.modelav.plot <- ggplot(data = dplyr::filter(time.modelav.summary.grouped,
                                           year >= 2016.75),
                      aes(x = year, y = model.averaged.relreduc)) +
  geom_line()
plot(impact.modelav.plot)


## Time-course of HIV incidence under the 2 scenarios (model-averaged)
inc.modelav.plot <- ggplot(data = input.output.F.long.revised.inc.average,
          aes(x = year,
              y = inc.av,
              group = scenario,
              col = scenario.av)) +
  geom_line(size = 1.5) +
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
                                "Counterfactual model 3"))
plot(inc.modelav.plot)

input.output.F.long.revised.inc.average[input.output.F.long.revised.inc.average$year == 2016.75, ]
input.output.F.long.revised.inc.average[input.output.F.long.revised.inc.average$year == 2020, ]


#### HIV PREVALENCE PLOT

# Overlaying UNAIDS HIV prevalence estimates
unaids.prev <- read.csv(file = "Eswatini_UNAIDS_2018_People living with HIV_National.csv",
                           sep = ";",
                           skip = 2,
                           header = FALSE,
                           col.names = c("year",
                                         "prev.est",
                                         "prev.lower",
                                         "prev.upper"))
unaids.art.cov$unique.id <- NA
unaids.art.cov$scenario <- NA


# 180 incidence estimates + # from = 7.25, to = 52, by = 0.25
# 180 incident cases +
# 129 art cases +
# 1 SHIMS1 prev +
# 1 SHIMS1 inc +
# 1 SHIMS2 inc +
# 92 annual HIV prev (from mid 1986 till 2032 every 6 months: 92 time points)

input.output.F.wide.revised.for.prev <- input.output.F.wide.revised[, c((x.offset+30+65+1+360+129+3):(x.offset+(30+65+360+129+3+92)), (index.end - 2):index.end)] # Making a copy for prevalence plot
names(input.output.F.wide.revised.for.prev)[1:92] <- seq(from = 1986.5,
                                                         to = 2032,
                                                         by = 0.5)

input.output.F.long.revised.prev <- gather(data = input.output.F.wide.revised.for.prev,
                                   key = year,
                                   value = hivprevalence,
                                   1:92)
input.output.F.long.revised.prev$year <- as.numeric(input.output.F.long.revised.prev$year)
input.output.F.long.revised.prev$scenario.model <- paste0(input.output.F.long.revised.prev$scenario,
                                                         input.output.F.long.revised.prev$model)

## Adding average HIV prevalence per time point and scenario, averaged over reps
input.output.F.long.revised.prev.average <- input.output.F.long.revised.prev %>%
  group_by(., scenario, year) %>%
  summarise(prev.av = mean(hivprevalence))
input.output.F.long.revised.prev.average$scenario.av <- paste0(input.output.F.long.revised.prev.average$scenario,
                                                               "overall average")
## Adding average HIV prevalence per time point, scenario, and model, averaged over reps
input.output.F.long.revised.prev.model.average <- input.output.F.long.revised.prev %>%
  group_by(., scenario, model, year) %>%
  summarise(prev.av = mean(hivprevalence))
input.output.F.long.revised.prev.model.average$scenario.model <- paste0(input.output.F.long.revised.prev.model.average$scenario,
                                                                       input.output.F.long.revised.prev.model.average$model)


# FaFc prevalence plot
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
  theme(panel.background = element_rect(fill = "grey97"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

plot(prev.FaFc.revised.plot)
# Saving the prevalence plot
ggsave(filename = "prev_FaFc_3models_revised_plot.pdf",
       plot = prev.FaFc.revised.plot,
       path = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/plots",
       width = 16, height = 10, units = "cm")
