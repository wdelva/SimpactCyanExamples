#### Run time analysis, as a function of initial number of individuals, events

library(RSimpactCyan)
library(RSimpactHelper)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(metafolio)


# Study design
# pop size: 1000, 2000, 5000, 10000,
# population.eyecap.fraction: 0.1, 0.2, 0.4,
# parallel: yes, no,
# concurrency: inputvector[12] and inputvector[13]: -0.5 and -5
# 10 sims per scenario
# 4 by 3 by 2 by 2 by 10 = 480

study.design.df <- as.data.frame(expand.grid(pop.size = c(2000, 5000, 10000),
                                             eyecap.fraction = c(0.2, 0.4),
                                             parallel.boolean = c(TRUE, FALSE),
                                             numrel.effect = c(-0.5, -5),
                                             seed.id = 1:10,
                                             KEEP.OUT.ATTRS = FALSE))
study.obs <- nrow(study.design.df)

study.output.df <- data.frame(sim.duration = rep(NA, study.obs),
                              sim.duration.units = rep(NA, study.obs),
                              sim.events = rep(NA, study.obs),
                              sim.rels = rep(NA, study.obs))

for (sim.index in 1:nrow(study.design.df)){
  
  
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
                                -0.5,
                                -2.5,
                                -1,
                                -5,
                                2,
                                2,
                                1,
                                0)
  
  inputvector <- c(study.design.df$seed.id[sim.index], inputvector.EAAA.example)
  inputvector[12] <- study.design.df$numrel.effect[sim.index]
  inputvector[13] <- study.design.df$numrel.effect[sim.index]
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  cfg.list <- input.params.creator(population.eyecap.fraction = study.design.df$eyecap.fraction[sim.index], # 0.2, 
                                   population.simtime = 40,
                                   population.nummen = study.design.df$pop.size[sim.index], # 2000,
                                   population.numwomen = study.design.df$pop.size[sim.index], # 2000, 
                                   population.msm = "no", hivseed.time = 10, hivseed.type = "amount", 
                                   hivseed.amount = 20, hivseed.age.min = 20, hivseed.age.max = 50, 
                                   hivtransmission.param.a = -1, hivtransmission.param.b = -90, 
                                   hivtransmission.param.c = 0.5, hivtransmission.param.f1 = log(2), 
                                   hivtransmission.param.f2 = log(log(sqrt(2))/log(2))/5, 
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01, 
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01, 
                                   formation.hazard.agegapry.meanage = -0.025, formation.hazard.agegapry.gap_factor_man_const = 0, 
                                   formation.hazard.agegapry.gap_factor_woman_const = 0, 
                                   formation.hazard.agegapry.gap_factor_man_exp = -1, formation.hazard.agegapry.gap_factor_woman_exp = -1, 
                                   formation.hazard.agegapry.gap_agescale_man = 0.25, formation.hazard.agegapry.gap_agescale_woman = 0.25, 
                                   dissolution.alpha_4 = -0.05, debut.debutage = 15, conception.alpha_base = -2.7, 
                                   dropout.interval.dist.type = "exponential")
  mu.cd4 <- 800
  var.cd4 <- 200^2
  mu.cd4.end <- 20
  var.cd4.end <- 5
  cfg.list["person.cd4.start.dist.type"] <- "lognormal"
  cfg.list["person.cd4.start.dist.lognormal.zeta"] <- log(mu.cd4/sqrt(1 + 
                                                                        var.cd4/mu.cd4^2))
  cfg.list["person.cd4.start.dist.lognormal.sigma"] <- sqrt(log(1 + 
                                                                  var.cd4/mu.cd4^2))
  cfg.list["person.cd4.end.dist.type"] <- "lognormal"
  cfg.list["person.cd4.end.dist.lognormal.zeta"] <- log(mu.cd4.end/sqrt(1 + 
                                                                          var.cd4.end/mu.cd4.end^2))
  cfg.list["person.cd4.end.dist.lognormal.sigma"] <- sqrt(log(1 + 
                                                                var.cd4.end/mu.cd4.end^2))
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0.3
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  cfg.list["person.agegap.man.dist.type"] <- "normal"
  cfg.list["person.agegap.woman.dist.type"] <- "normal"
  cfg.list["monitoring.cd4.threshold"] <- 1
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
  cfg.list["diagnosis.baseline"] <- -99999
  cfg.list["periodiclogging.interval"] <- 0.25
  cfg.list["dropout.interval.dist.exponential.lambda"] <- 0.1
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * 
    as.numeric(cfg.list["population.nummen"][1]) * 6
  cfg.list["person.vsp.toacute.x"] <- 5
  seedid <- inputvector[1]
  cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
  cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2]))/log(inputvector[2]))/5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]
  cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10]
  cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10]
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[11]
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
  cfg.list["conception.alpha_base"] <- inputvector[14]
  cfg.list["dissolution.alpha_0"] <- inputvector[15]
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- inputvector[16]
  art.intro["monitoring.cd4.threshold"] <- 100
  art.intro["formation.hazard.agegapry.baseline"] <- inputvector[11] - 
    0.5
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- inputvector[16] + inputvector[17]
  art.intro1["monitoring.cd4.threshold"] <- 150
  art.intro2 <- list()
  art.intro2["time"] <- 23
  art.intro2["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + 
    inputvector[18]
  art.intro2["monitoring.cd4.threshold"] <- 200
  art.intro2["formation.hazard.agegapry.baseline"] <- inputvector[11] - 
    1
  art.intro3 <- list()
  art.intro3["time"] <- 30
  art.intro3["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + 
    inputvector[18] + inputvector[19]
  art.intro3["monitoring.cd4.threshold"] <- 350
  art.intro4 <- list()
  art.intro4["time"] <- 33.5
  art.intro4["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + 
    inputvector[18] + inputvector[19] + inputvector[20]
  art.intro4["monitoring.cd4.threshold"] <- 500
  art.intro5 <- list()
  art.intro5["time"] <- 36.75
  art.intro5["monitoring.cd4.threshold"] <- 6000
  ART.factual <- list(art.intro, art.intro1, art.intro2, art.intro3, 
                      art.intro4, art.intro5)
  ART.counterfactual <- list(art.intro, art.intro1, art.intro2, 
                             art.intro3)
  identifier <- paste0(seedid)
  rootDir <- "/tmp"
  destDir <- paste0(rootDir, "/", identifier)
  
  
  
  sim.start <- Sys.time()
  results <- tryCatch(simpact.run(configParams = cfg.list, 
                                  destDir = destDir, agedist = age.distr, intervention = ART.factual, 
                                  seed = seedid, identifierFormat = identifier,
                                  parallel = study.design.df$parallel.boolean[sim.index]), 
                      error = simpact.errFunction)
  sim.end <- Sys.time()
  
  # Duration of simulation
  sim.duration <- sim.end - sim.start
  
  # Number of events
  sim.events <- as.numeric(results["eventsexecuted"])
  
  if (length(results) > 0){
    # Number of relationships
    datalist.EAAA <- readthedata(results)
    agemix.episodes.df <- agemix.episodes.df.maker(datalist.EAAA)
    agemix.rels.df <- agemix.rels.df.maker(dataframe = agemix.episodes.df,
                                           agegroup = c(15, 200),
                                           timewindow = 40,
                                           timepoint = 40, start = FALSE)
    sim.rels <- as.numeric(nrow(agemix.rels.df))
    
    study.output.df$sim.duration[sim.index] <- sim.duration
    study.output.df$sim.duration.units[sim.index] <- units(sim.duration)
    study.output.df$sim.events[sim.index] <- sim.events
    study.output.df$sim.rels[sim.index] <- sim.rels
  }

  
  # cleaning up
  unlink(paste0(rootDir, "/", identifier), recursive = TRUE)
}

study.data <- cbind(study.design.df, study.output.df)
save(study.data, file = "/Users/wimdelva/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/study.data.RData")
table(study.data$sim.duration.units)
study.data$sim.duration.secs <- study.data$sim.duration
study.data$sim.duration.secs[study.data$sim.duration.units == "mins"] <- study.data$sim.duration[study.data$sim.duration.units == "mins"] * 60

study.data.grouped <- group_by(study.data, pop.size, eyecap.fraction, numrel.effect, parallel.boolean)
mean.rels.grouped <- summarize(study.data.grouped,
                               mean.rels = mean(sim.rels, na.rm = T),
                               mean.events = mean(sim.events, na.rm = T),
                               mean.runtime = mean(sim.duration.secs, na.rm = T))

parallel.labels <- c("TRUE" = "Multi-core", "FALSE" = "Single-core")

runtime.plot <- ggplot(data = mean.rels.grouped,
                       aes(x = 2 * pop.size,
                           y = mean.runtime / 60,
                           col = mean.rels,
                           shape = factor(eyecap.fraction))) +
  geom_point(size = 3) +
  facet_wrap(~parallel.boolean, labeller = labeller(parallel.boolean = parallel.labels)) +
  scale_colour_distiller("Mean number of\nrelationships formed",
                         palette = "Spectral") +
  scale_shape_discrete(name = "Accessible\npopulation fraction") +
  xlab("Population size") +
  ylab("Mean runtime (minutes)") +
  theme(panel.background = element_rect(fill = "grey97"))
runtime.plot
