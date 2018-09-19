setwd("/path/to/your/working_directory/") # Change this to your own working directory

library(RSimpactCyan)
library(RSimpactHelper)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(metafolio)

inputvector.EAAA.example <- c(1.1,
                                 0.4,
                                 1,
                                 3.5,
                                 0.8,
                                 0.3,
                                 45,
                                 30,
                                 -1,
                                 3.8,
                                 -0.5,
                                 -0.45,
                                 -2.5,
                                 -1,
                                 -5,
                                 2,
                                 2,
                                 1,
                                 0)

reps <- 20

EAAA.design.matrix <- matrix(rep(inputvector.EAAA.example,
                                 times = reps),
                             nrow = reps,
                             byrow = TRUE)

EAAA.revised.F.a.output <- simpact.parallel(model = EAAA.revised.F.a.wrapper,
                                            actual.input.matrix = EAAA.design.matrix,
                                            seed_count = 0,
                                            n_cluster = 4)

EAAA.revised.F.c.output <- simpact.parallel(model = EAAA.revised.F.c.wrapper,
                                            actual.input.matrix = EAAA.design.matrix,
                                            seed_count = 0,
                                            n_cluster = 4)

input.df <- data.frame(seed = 1:reps,
                        EAAA.design.matrix)

input.output.Fa.df <- cbind(input.df, as.data.frame(EAAA.revised.F.a.output))
input.output.Fc.df <- cbind(input.df, as.data.frame(EAAA.revised.F.c.output))
input.output.Fa.df$scenario <- "Fa"
input.output.Fc.df$scenario <- "Fc"


## combining Fa and Fc datasets
input.output.F.wide <- rbind(input.output.Fa.df,
                             input.output.Fc.df)
input.output.F.wide$unique.id <- 1:nrow(input.output.F.wide)

save(input.output.F.wide, file = "/path/to/your/working_directory/SimpactPaperEAAAExample.RData")

x.offset <- 20 # seed + 19 parameters

# Colour scheme
# dark red #7D0E0F
# light red #E41A1C
# dark blue #20496B
# light blue #377EB8

cols <- c("#E41A1C", "#377EB8")
darkcols <- c("#7D0E0F", "#20496B")

#### ART COVERAGE PLOT
input.output.F.wide.for.cov <- input.output.F.wide[, c((x.offset+30):(x.offset+72), 616:617)] # Making a copy for ART coverage plot
names(input.output.F.wide.for.cov)[1:43] <- 1990:2032

input.output.F.long.cov <- gather(data = input.output.F.wide.for.cov,
                                  key = year,
                                  value = art.cov,
                                  1:43)
input.output.F.long.cov$year <- as.numeric(input.output.F.long.cov$year)

## Adding average ART coverage per year and scenario, averaged over reps
input.output.F.long.cov.average <- input.output.F.long.cov %>%
  group_by(., scenario, year) %>%
  summarise(cov.av = mean(art.cov))
input.output.F.long.cov.average$scenario.av <- paste0(input.output.F.long.cov.average$scenario, ".av")



# FaFc ART coverage plot
cov.FaFc.plot <- ggplot(data = input.output.F.long.cov,
                        aes(year, art.cov, group = unique.id, colour = scenario)) +
  geom_line(size = 0.3) +
  theme(legend.position=c(0.72, 0.52),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  xlab("Time (years)") +
  ylab("ART coverage") +
  ylim(0, 1) +
  geom_line(data = input.output.F.long.cov.average,
             aes(x = year, y = cov.av, group = scenario, col = scenario.av),
             size = 1) +
  scale_color_manual(values = c("Fa" = cols[1],
                                "Fc" = cols[2],
                                "Fa.av" = darkcols[1],
                                "Fc.av" = darkcols[2]),
                     guide = FALSE)
plot(cov.FaFc.plot)
# Saving the ART coverage plot
ggsave(filename = "cov.FaFc.plot.pdf",
       plot = cov.FaFc.plot,
       path = "/path/to/your/working_directory/plots",
       width = 16, height = 10, units = "cm")


#### HIV INCIDENCE PLOT
input.output.F.wide.for.inc <- input.output.F.wide[, c((x.offset+73+1):(x.offset+73+168), 616:617)] # Making a copy for incidence plot
names(input.output.F.wide.for.inc)[1:168] <- seq(1990, 2031.75, by = 0.25)

input.output.F.long.inc <- gather(data = input.output.F.wide.for.inc,
                                  key = year,
                                  value = hivincidence,
                                  1:168)
input.output.F.long.inc$year <- as.numeric(input.output.F.long.inc$year)

## Adding average HIV incidence per time point and scenario, averaged over reps
input.output.F.long.inc.average <- input.output.F.long.inc %>%
  group_by(., scenario, year) %>%
  summarise(inc.av = mean(hivincidence))
input.output.F.long.inc.average$scenario.av <- paste0(input.output.F.long.inc.average$scenario, ".av")


# FaFc incidence plot
inc.FaFc.plot <- ggplot(data = input.output.F.long.inc,
                        aes(year, hivincidence, group = unique.id, colour = scenario)) +
  geom_line(size = 0.2) +
  xlab("Time (years)") +
  ylab("HIV incidence") +
  geom_line(data = input.output.F.long.inc.average,
            aes(x = year, y = inc.av, group = scenario, col = scenario.av),
            size = 1) +
  scale_color_manual(values = c("Fa" = cols[1],
                                "Fc" = cols[2],
                                "Fa.av" = darkcols[1],
                                "Fc.av" = darkcols[2]),
                     guide = FALSE)
plot(inc.FaFc.plot)
# Saving the incidence plot
ggsave(filename = "inc.FaFc.plot.pdf",
       plot = inc.FaFc.plot,
       path = "/path/to/your/working_directory/plots",
       width = 16, height = 10, units = "cm")


#### HIV PREVALENCE PLOT
input.output.F.wide.for.prev <- input.output.F.wide[, c(562:614, 616:617)] # Making a copy for prevalence plot
names(input.output.F.wide.for.prev)[1:53] <- seq(1980, 2032, by = 1)

input.output.F.long.prev <- gather(data = input.output.F.wide.for.prev,
                                  key = year,
                                  value = hivprevalence,
                                  1:53)
input.output.F.long.prev$year <- as.numeric(input.output.F.long.prev$year)

## Adding average HIV prevalence per time point and scenario, averaged over reps
input.output.F.long.prev.average <- input.output.F.long.prev %>%
  group_by(., scenario, year) %>%
  summarise(prev.av = mean(hivprevalence))
input.output.F.long.prev.average$scenario.av <- paste0(input.output.F.long.prev.average$scenario, ".av")


# FaFc prevalence plot
prev.FaFc.plot <- ggplot(data = input.output.F.long.prev,
                        aes(year, hivprevalence, group = unique.id, colour = scenario)) +
  geom_line(size = 0.2) +
  xlab("Time (years)") +
  ylab("HIV prevalence") +
  geom_line(data = input.output.F.long.prev.average,
            aes(x = year, y = prev.av, group = scenario, col = scenario.av),
            size = 1) +
  scale_color_manual(values = c("Fa" = cols[1],
                                "Fc" = cols[2],
                                "Fa.av" = darkcols[1],
                                "Fc.av" = darkcols[2]),
                     guide = FALSE)
plot(prev.FaFc.plot)
# Saving the prevalence plot
ggsave(filename = "prev.FaFc.plot.pdf",
       plot = prev.FaFc.plot,
       path = "/path/to/your/working_directory/plots",
       width = 16, height = 10, units = "cm")

# HIV prevalence among 18-50 year olds at the time of the SHIMS I survey:
input.output.F.wide[, c(559, 616)] %>%
  group_by(., scenario) %>%
  summarise(prev.av = mean(V539))


