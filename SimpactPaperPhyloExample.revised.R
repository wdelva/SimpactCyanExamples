# Before running the code below, make sure you have followed the instructions for installing all the required software,
# as explained in the README file.
setwd("/Users/wimdelva/Documents/SimpactCyanExamples")

# Make sure compiled tools (Seq-Gen and FastTree) are in same working directory

# Enable installing from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("ggtree")
## Load required packages
library(RSimpactCyan)
library(devtools)
install_github("wdelva/RSimpactHelp")
library(RSimpactHelper)
library(tidyverse)
library(ggplot2)

# Sourcing the IDs.Seq.Random.skew function and renaming IDS in transmission table 
# to match tips names in the tree
source("util.seq.cov.weight.R")

#######################
# Step 1: Run Simpact #
#######################

EAAA.SciRep.revision.Seq.loaded.object <- load(file = paste0(getwd(), "/EAAA.SciRep.revision2.Seq.RData"))
# load(file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/EAAA.SciRep.revision2.Seq.RData")

Seq.object <- EAAA.SciRep.revision2.Seq 
# The vector of target features was:
features.pop.growth <- exp(0.015)
features.hiv.prev <- c(0.143, 0.008, 0.315, 0.066, 0.467, 0.213, 0.538, 0.366, 0.491, 0.47, 0.397, 0.455, 0.316, 0.425)
features.hiv.inc <- exp(c(0.038, 0.008, 0.043, 0.016, 0.02, 0.026, 0.027, 0.031, 0.04, 0.004, 0.021, 0.012, 0.012, 0))
features.art.cov <- c(0.37, 0.40, 0.44, 0.49, 0.58, 0.67, 0.76, 0.85) # c(0.33, 0.38, 0.45, 0.51, 0.61, 0.7, 0.8)
features.vl.suppr <- 0.74  # 0.68
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


# Using the best-fitting model from the calibration
inputvector.EAAA.example <- Seq.object$param[index.bestfit, ]


inputvector <- c(0, inputvector.EAAA.example)

age.distr <- agedistr.creator(shape = 5, scale = 65)

cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 52,
                                 population.nummen = 5000,
                                 population.numwomen = 5000,
                                 population.msm = "no",
                                 hivseed.time = 8.5,
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 hivtransmission.param.a = -1,
                                 hivtransmission.param.b = -90,
                                 hivtransmission.param.c = 0.5,
                                 hivtransmission.param.f1 = log(2),
                                 hivtransmission.param.f2 = log(log(sqrt(2)) / log(2)) / 5,
                                 formation.hazard.agegapry.gap_factor_man_age = -0.01,
                                 formation.hazard.agegapry.gap_factor_woman_age = -0.01,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_factor_man_exp = -1,
                                 formation.hazard.agegapry.gap_factor_woman_exp = -1,
                                 formation.hazard.agegapry.gap_agescale_man = 0.25,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.25,
                                 dissolution.alpha_4 = -0.05,
                                 debut.debutage = 15,
                                 conception.alpha_base = -2.7,
                                 dropout.interval.dist.type = "uniform")

#standard deviation of 200 CD4 cells
#mu = ln(mean / sqrt(1 + variance/mean^2))
#sigma^2 = ln(1 + variance/mean^2)
#Here, we say mean = 825 and variance = 200^2
mu.cd4 <- 800
var.cd4 <- 200^2
mu.cd4.end <- 20
var.cd4.end <- 5
cfg.list["person.cd4.start.dist.type"] <- "lognormal"
cfg.list["person.cd4.start.dist.lognormal.zeta"] <- log(mu.cd4/sqrt(1+var.cd4/mu.cd4^2))
cfg.list["person.cd4.start.dist.lognormal.sigma"] <- sqrt(log(1+var.cd4/mu.cd4^2))
cfg.list["person.cd4.end.dist.type"] <- "lognormal"
cfg.list["person.cd4.end.dist.lognormal.zeta"] <- log(mu.cd4.end/sqrt(1+var.cd4.end/mu.cd4.end^2))
cfg.list["person.cd4.end.dist.lognormal.sigma"] <- sqrt(log(1+var.cd4.end/mu.cd4.end^2))

cfg.list["formation.hazard.agegapry.baseline"] <- 2
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3

cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal"
cfg.list["person.agegap.woman.dist.type"] <- "normal"

cfg.list["monitoring.cd4.threshold"] <- 1 # 0
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75 # 1 # 0.9 # 0.75 # 0.5
cfg.list["diagnosis.baseline"] <- -99999 # -2
cfg.list["periodiclogging.interval"] <- 0.25
# cfg.list["dropout.interval.dist.exponential.lambda"] <- 0.1


cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 6

cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1]
cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
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

# Introducing ART
art.intro <- list()
art.intro["time"] <- 20
art.intro["diagnosis.baseline"] <- inputvector[16] # prior [-4 , 0] # -2
art.intro["monitoring.cd4.threshold"] <- 100

art.intro1 <- list()
art.intro1["time"] <- 22
art.intro1["diagnosis.baseline"] <- inputvector[16] + inputvector[17] # prior [0, 2] # -1.8
art.intro1["monitoring.cd4.threshold"] <- 150

art.intro2 <- list()
art.intro2["time"] <- 23
art.intro2["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + inputvector[18] # prior [0, 2] # -1.5
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 30
art.intro3["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + inputvector[18] + inputvector[19] # prior [0, 2] #-1
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 33.5
art.intro4["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + inputvector[18] + inputvector[19] + inputvector[20] # prior [0, 2]
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 36.75
art.intro5["monitoring.cd4.threshold"] <- 6000


ART.factual <- list(art.intro,art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)

identifier <- paste0(seedid)
destDir <- paste0(getwd(), "/temp")

results <- tryCatch(simpact.run(configParams = cfg.list,
                                destDir = destDir,
                                agedist = age.distr,
                                intervention = ART.factual,
                                seed = seedid),
                    error = simpact.errFunction)
datalist <- readthedata(results)

###########################################
# Step 2: Construct transmission networks #
###########################################
# Produce a list of transmission networks in epi object format
simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)

net.sizes <- purrr::map(simpact.trans.net, 1) %>%
  lapply(., length) %>%
  unlist()
max.net.size <- max(net.sizes)
max.net.size.index <- which(net.sizes %in% max.net.size) # So we know which network and tree to visualise

###############################
# Step 3: Sequence simulation #
###############################

# To do the forward simulation, we need to specify the molecular evolution model.
# To inform this, we did a phylogenetic analysis of the polCDS gene region of 386 sequences from patients with HIV-1 subtype C infection from South Africa (one sequence per patient). Data retrieved from:
# https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html

# These sequences were analysed with jModelTest version 2.1.3, downloadable at: https://github.com/ddarriba/jmodeltest2

hiv_db.df <- phylotools::read.fasta(file = "hiv_db.fasta")
hiv_db.df$seq.name <- gsub('^..', '', hiv_db.df$seq.name)

hiv_za.df <- hiv_db.df[grep("^ZA", hiv_db.df$seq.name), ]
phylotools::dat2fasta(hiv_za.df, outfile = "hiv_za.fasta")

dirseqgen <- getwd()
dirfasttree <- getwd()
sub.dir.rename <- getwd()

sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                               sub.dir.rename = sub.dir.rename,
                               simpact.trans.net = simpact.trans.net,
                               seq.gen.tool = "seq-gen",
                               seeds.num = 777,
                               endpoint = 40,
                               limitTransmEvents = max.net.size,
                               hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                               clust = FALSE)

# Output sequence file: C.Epidemic_seed.seq.bis.sim.nwk.fasta

# Transform the sequence format to be handled by ClusterPicker
sequ.dna <- ape::read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
ape::write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")


#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################

fitGTRGI.top <- phylogenetic.tree.phangorn.par(simseqfile = sequ.dna)
unrooted.tree <- fitGTRGI.top$tree

# Before we calibrate this tree, we need to root it.
calendar.dates = paste0("samplingtimes_seed_", max.net.size.index, ".csv")
samp.dates <- read.csv(paste0(sub.dir.rename, "/", calendar.dates))
time.samp <- dates.Transform.NamedVector(dates = samp.dates,
                                         count.start = 1977,
                                         endsim = 40) # name the dates


### Match dates and phylogenetic tree leaves ###
################################################
time.samp.df <- data.frame(samp.ID = names(time.samp),
                           time.samp = time.samp)
tree.tips.df <- data.frame(samp.ID = unrooted.tree$tip.label)
Ord.tree.dates <- left_join(x = tree.tips.df,
          y = time.samp.df) %>%
           dplyr::select(time.samp) %>%
           unlist()
names(Ord.tree.dates) <- tree.tips.df$samp.ID

rooted.tree <- ape::rtt(t = unrooted.tree,
         tip.dates = Ord.tree.dates,
         ncpu = 1,
         objective = "correlation")


# Calibrate the phylogenetic tree
dater.tree <- treedater::dater(rooted.tree, Ord.tree.dates, s = 3000, searchRoot = 100) # s is the length of sequence


# Node age with picante package
N <- picante::node.age(dater.tree)

# Time to MRCA: internal nodes ages
int.node.age <- N$Ti

latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date

#####################################################################################
# Step 5: Compute transmission events and internal nodes in one year time intervals #
#####################################################################################

## Entire transmission network 
trans.net <- tips.labels(simpact.trans.net = simpact.trans.net,
                                     limitTransmEvents = max.net.size)

trans.net$dtimes <- abs(trans.net$dtimes-40)+1977
trans.net$itimes <- abs(trans.net$itimes-40)+1977
trans.net$DonId <- factor(trans.net$DonId)
trans.net$RecId <- factor(trans.net$RecId)


gag <- as.matrix(trans.net)
# gag <- gag[-1,] # remove universall seed -1
ga.graph = igraph::graph.edgelist(gag[ , 7:6])

edges <- as.data.frame(igraph::get.edgelist(ga.graph))
names(edges) <- c("infector", "infectee")
vertices <- as.data.frame(igraph::vertex_attr(ga.graph))
network.df <- fortify(geomnet::as.edgedf(edges), vertices)

#########################################################################
# Step 6: Creating the panels A, B and C of the figure for this example #
#########################################################################

# A. Transmission network
# The transmission event from the "ghost infector with ID "-1" to the seed individual with ID "1526" is excluded from the transmission network
transmissionnetwork.plot <- ggplot(data = network.df[-1, ],
                                   aes(from_id = from_id,
                                       to_id = to_id)) +
  geomnet::geom_net(directed = TRUE,
           size = 1, # 2.5
           layout.alg = "kamadakawai",
           alpha = 1,
           arrowsize = 0.15, # 0.5
           arrowgap = 0.01,
           ecolour = "darkgrey",
           colour = "black") +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  ylab("")
print(transmissionnetwork.plot)

ggsave(filename = "network.pdf",
       plot = transmissionnetwork.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")



# B. Phylogenetic tree

tree <- dater.tree
class(tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
mrsd <- max(dater.tree$sts)

dates <- format(lubridate::date_decimal(mrsd), "%Y-%m-%d")


# Adding the root edge
# Seed individual was introduced in 1985.5
root.edge.length <- min(N$Ti) - 1985.5
tree <- TreePar::addroot(tree, root.edge.length)

phylotree.plot <- ggtree::ggtree(tree,
                         mrsd = dates[1],
                         size = 0.05) + 
  ggtree::theme_tree2() +
  theme_grey() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) +
  xlab("Time") +
  ylab("")
print(phylotree.plot)

ggsave(filename = "tree.pdf",
       plot = phylotree.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")




# C. Transmission event versus internal nodes
# The transmission event from the "ghost infector with ID "-1" to the seed individual with ID "1526" is excluded from the transmission network
annual.infections <- as.numeric(table(floor(trans.net$itimes[-1])))
year.bins <- as.numeric(names(table(floor(trans.net$itimes[-1]))))
transmission.events.df <- data.frame(number.trans = annual.infections,
                                     calendaryear = year.bins)
                                     
# Timing of internal nodes
aged.dater.tree <- picante::node.age(dater.tree)
annual.internal.nodes <- as.numeric(table(floor(aged.dater.tree$Ti)))
year.bins.intnodes <- as.numeric(names(table(floor(aged.dater.tree$Ti))))
internal.nodes.df <- data.frame(number.nodes = annual.internal.nodes,
                                     calendaryear = year.bins.intnodes)

years.df <- data.frame(calendaryear = 1985:2017)
trans.df <- dplyr::left_join(x = years.df,
                 y = transmission.events.df) %>%
  replace_na(list(number.trans = 0))
trans.and.nodes.df <- dplyr::left_join(x = trans.df,
                                       y = internal.nodes.df) %>%
  replace_na(list(number.nodes = 0))


trans.and.nodes.long.df <- gather(trans.and.nodes.df,
                                  key = "Events",
                                  value = "Number",
                                  number.trans:number.nodes,
                                  factor_key = TRUE)


transandnodes.plot <- ggplot(data = dplyr::filter(trans.and.nodes.long.df,
                                                  calendaryear >= 1985,
                                                  calendaryear <= 2017),
                             aes(x = calendaryear,
                                 y = Number,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time")
print(transandnodes.plot)

ggsave(filename = "transandnodes.pdf",
       plot = transandnodes.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")

## Adding density to dataset trans.and.nodes.long.df
trans.and.nodes.long.enriched.df <- trans.and.nodes.long.df %>%
  group_by(Events) %>%
  mutate(total.events = sum(Number),
         percentage.events = Number / total.events)

# New plot comparing densities  
transandnodesfraction.plot <- ggplot(data = dplyr::filter(trans.and.nodes.long.enriched.df,
                                     calendaryear >= 1985,
                                     calendaryear <= 2017),
                             aes(x = calendaryear,
                                 y = percentage.events,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time") +
  ylab("Fraction")
print(transandnodesfraction.plot)

ggsave(filename = "transandnodesfraction.pdf",
       plot = transandnodesfraction.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")




#########################################
##### Sequence coverage scenarios #######
#########################################


# Scenario 1: where we create a dataset with only 50% of the sequences, 
# and they are sampled completely at random. ------------------------------


# Select IDs 
#############

seq.cov <- 0.5
set.seed(0)
sampled.50.random.indices <- sample.int(nrow(sequ.dna), size = round(nrow(sequ.dna) * seq.cov), replace = FALSE)
sequ.dna.50.random <- sequ.dna[sampled.50.random.indices, ]


# Build and calibrate the phylogenetic tree
############################################

fitGTRGI.top.50.random <- phylogenetic.tree.phangorn.par(simseqfile = sequ.dna.50.random)
unrooted.tree.50.random <- fitGTRGI.top.50.random$tree


### Match dates and phylogenetic tree leaves ###
################################################
time.samp.df <- data.frame(samp.ID = names(time.samp),
                           time.samp = time.samp)
tree.tips.50.random.df <- data.frame(samp.ID = unrooted.tree.50.random$tip.label)
Ord.tree.dates.50.random <- left_join(x = tree.tips.50.random.df,
                            y = time.samp.df) %>%
  dplyr::select(time.samp) %>%
  unlist()
names(Ord.tree.dates.50.random) <- tree.tips.50.random.df$samp.ID


rooted.tree.50.random <- ape::rtt(t = unrooted.tree.50.random,
                        tip.dates = Ord.tree.dates.50.random,
                        ncpu = 1,
                        objective = "correlation")


# Calibrate the phylogenetic tree
dater.tree.50.random <- treedater::dater(rooted.tree.50.random, Ord.tree.dates.50.random, s = 3000, searchRoot = 100) # s is the length of sequence

# write.tree(tree.calib, file = paste0(sub.dir.rename, paste0("/calibrated.tree.",simseqfile,".tree")))
# calibrated.tree.C.Epidemic.fas.tree

# Node age with picante package
N.50.random <- picante::node.age(dater.tree.50.random)

# Time to MRCA: internal nodes ages
int.node.age.50.random <- N.50.random$Ti

latest.samp.50.random <- N.50.random$timeToMRCA + N.50.random$timeOfMRCA # latest sampling date


# C. Transmission event versus internal nodes
# The transmission event from the "ghost infector with ID "-1" to the seed individual with ID "1526" is excluded from the transmission network
annual.infections <- as.numeric(table(floor(trans.net$itimes[-1])))
year.bins <- as.numeric(names(table(floor(trans.net$itimes[-1]))))
transmission.events.df <- data.frame(number.trans = annual.infections,
                                     calendaryear = year.bins)

# Timing of internal nodes
aged.dater.tree.50.random <- picante::node.age(dater.tree.50.random)
annual.internal.nodes.50.random <- as.numeric(table(floor(aged.dater.tree.50.random$Ti)))
year.bins.intnodes.50.random <- as.numeric(names(table(floor(aged.dater.tree.50.random$Ti))))
internal.nodes.50.random.df <- data.frame(number.nodes = annual.internal.nodes.50.random,
                                calendaryear = year.bins.intnodes.50.random)

years.df <- data.frame(calendaryear = 1985:2017)
trans.df <- dplyr::left_join(x = years.df,
                             y = transmission.events.df) %>%
  replace_na(list(number.trans = 0))
trans.and.nodes.50.random.df <- dplyr::left_join(x = trans.df,
                                       y = internal.nodes.50.random.df) %>%
  replace_na(list(number.nodes = 0))


trans.and.nodes.long.50.random.df <- gather(trans.and.nodes.50.random.df,
                                  key = "Events",
                                  value = "Number",
                                  number.trans:number.nodes,
                                  factor_key = TRUE)


transandnodes.50.random.plot <- ggplot(data = trans.and.nodes.long.50.random.df,
                             aes(x = calendaryear,
                                 y = Number,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time")
print(transandnodes.50.random.plot)

ggsave(filename = "transandnodes.50.random.pdf",
       plot = transandnodes.50.random.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")

## Adding density to dataset trans.and.nodes.long.df
trans.and.nodes.long.50.random.enriched.df <- trans.and.nodes.long.50.random.df %>%
  group_by(Events) %>%
  mutate(total.events = sum(Number),
         percentage.events = Number / total.events)

# New plot comparing densities  
transandnodes.fraction.50.random.plot <- ggplot(data = dplyr::filter(trans.and.nodes.long.50.random.enriched.df,
                                                                    calendaryear >= 1985,
                                                                    calendaryear <= 2017),
                             aes(x = calendaryear,
                                 y = percentage.events,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time") +
  ylab("Fraction")
print(transandnodes.fraction.50.random.plot)

ggsave(filename = "transandnodesfraction50random.pdf",
       plot = transandnodes.fraction.50.random.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")


####################################################
## Scenario 2: sample 25% of the sequences
####################################################

# Select IDs 
#############

seq.cov <- 0.25
set.seed(0)
sampled.25.random.indices <- sample.int(nrow(sequ.dna), size = round(nrow(sequ.dna) * seq.cov), replace = FALSE)
sequ.dna.25.random <- sequ.dna[sampled.25.random.indices, ]


# Build and calibrate the phylogenetic tree
############################################

fitGTRGI.top.25.random <- phylogenetic.tree.phangorn.par(simseqfile = sequ.dna.25.random)
unrooted.tree.25.random <- fitGTRGI.top.25.random$tree


### Match dates and phylogenetic tree leaves ###
################################################
time.samp.df <- data.frame(samp.ID = names(time.samp),
                           time.samp = time.samp)
tree.tips.25.random.df <- data.frame(samp.ID = unrooted.tree.25.random$tip.label)
Ord.tree.dates.25.random <- left_join(x = tree.tips.25.random.df,
                                      y = time.samp.df) %>%
  dplyr::select(time.samp) %>%
  unlist()
names(Ord.tree.dates.25.random) <- tree.tips.25.random.df$samp.ID


rooted.tree.25.random <- ape::rtt(t = unrooted.tree.25.random,
                                  tip.dates = Ord.tree.dates.25.random,
                                  ncpu = 1,
                                  objective = "correlation")


# Calibrate the phylogenetic tree
dater.tree.25.random <- treedater::dater(rooted.tree.25.random, Ord.tree.dates.25.random, s = 3000, searchRoot = 100) # s is the length of sequence

# write.tree(tree.calib, file = paste0(sub.dir.rename, paste0("/calibrated.tree.",simseqfile,".tree")))
# calibrated.tree.C.Epidemic.fas.tree

# Node age with picante package
N.25.random <- picante::node.age(dater.tree.25.random)

# Time to MRCA: internal nodes ages
int.node.age.25.random <- N.25.random$Ti

latest.samp.25.random <- N.25.random$timeToMRCA + N.25.random$timeOfMRCA # latest sampling date


# C. Transmission event versus internal nodes
# The transmission event from the "ghost infector with ID "-1" to the seed individual with ID "1526" is excluded from the transmission network
annual.infections <- as.numeric(table(floor(trans.net$itimes[-1])))
year.bins <- as.numeric(names(table(floor(trans.net$itimes[-1]))))
transmission.events.df <- data.frame(number.trans = annual.infections,
                                     calendaryear = year.bins)

# Timing of internal nodes
aged.dater.tree.25.random <- picante::node.age(dater.tree.25.random)
annual.internal.nodes.25.random <- as.numeric(table(floor(aged.dater.tree.25.random$Ti)))
year.bins.intnodes.25.random <- as.numeric(names(table(floor(aged.dater.tree.25.random$Ti))))
internal.nodes.25.random.df <- data.frame(number.nodes = annual.internal.nodes.25.random,
                                          calendaryear = year.bins.intnodes.25.random)

years.df <- data.frame(calendaryear = 1985:2017)
trans.df <- dplyr::left_join(x = years.df,
                             y = transmission.events.df) %>%
  replace_na(list(number.trans = 0))
trans.and.nodes.25.random.df <- dplyr::left_join(x = trans.df,
                                                 y = internal.nodes.25.random.df) %>%
  replace_na(list(number.nodes = 0))


trans.and.nodes.long.25.random.df <- gather(trans.and.nodes.25.random.df,
                                            key = "Events",
                                            value = "Number",
                                            number.trans:number.nodes,
                                            factor_key = TRUE)


transandnodes.25.random.plot <- ggplot(data = trans.and.nodes.long.25.random.df,
                                       aes(x = calendaryear,
                                           y = Number,
                                           colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time")
print(transandnodes.25.random.plot)

ggsave(filename = "transandnodes.25.random.pdf",
       plot = transandnodes.25.random.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")

## Adding density to dataset trans.and.nodes.long.df
trans.and.nodes.long.25.random.enriched.df <- trans.and.nodes.long.25.random.df %>%
  group_by(Events) %>%
  mutate(total.events = sum(Number),
         percentage.events = Number / total.events)

# New plot comparing densities  
transandnodes.fraction.25.random.plot <- ggplot(data = dplyr::filter(trans.and.nodes.long.25.random.enriched.df,
                                                                     calendaryear >= 1985,
                                                                     calendaryear <= 2017),
                                                aes(x = calendaryear,
                                                    y = percentage.events,
                                                    colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time") +
  ylab("Fraction")
print(transandnodes.fraction.25.random.plot)

ggsave(filename = "transandnodesfraction25random.pdf",
       plot = transandnodes.fraction.25.random.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")



####################################################
## Scenario 3: sample 12.5% of the sequences
####################################################

# Select IDs 
#############

seq.cov <- 0.125
set.seed(0)
sampled.125.random.indices <- sample.int(nrow(sequ.dna), size = round(nrow(sequ.dna) * seq.cov), replace = FALSE)
sequ.dna.125.random <- sequ.dna[sampled.125.random.indices, ]


# Build and calibrate the phylogenetic tree
############################################

fitGTRGI.top.125.random <- phylogenetic.tree.phangorn.par(simseqfile = sequ.dna.125.random)
unrooted.tree.125.random <- fitGTRGI.top.125.random$tree


### Match dates and phylogenetic tree leaves ###
################################################
time.samp.df <- data.frame(samp.ID = names(time.samp),
                           time.samp = time.samp)
tree.tips.125.random.df <- data.frame(samp.ID = unrooted.tree.125.random$tip.label)
Ord.tree.dates.125.random <- left_join(x = tree.tips.125.random.df,
                                      y = time.samp.df) %>%
  dplyr::select(time.samp) %>%
  unlist()
names(Ord.tree.dates.125.random) <- tree.tips.125.random.df$samp.ID


rooted.tree.125.random <- ape::rtt(t = unrooted.tree.125.random,
                                  tip.dates = Ord.tree.dates.125.random,
                                  ncpu = 1,
                                  objective = "correlation")


# Calibrate the phylogenetic tree
dater.tree.125.random <- treedater::dater(rooted.tree.125.random, Ord.tree.dates.125.random, s = 3000, searchRoot = 100) # s is the length of sequence

# write.tree(tree.calib, file = paste0(sub.dir.rename, paste0("/calibrated.tree.",simseqfile,".tree")))
# calibrated.tree.C.Epidemic.fas.tree

# Node age with picante package
N.125.random <- picante::node.age(dater.tree.125.random)

# Time to MRCA: internal nodes ages
int.node.age.125.random <- N.125.random$Ti

latest.samp.125.random <- N.125.random$timeToMRCA + N.125.random$timeOfMRCA # latest sampling date


# C. Transmission event versus internal nodes
# The transmission event from the "ghost infector with ID "-1" to the seed individual with ID "1526" is excluded from the transmission network
annual.infections <- as.numeric(table(floor(trans.net$itimes[-1])))
year.bins <- as.numeric(names(table(floor(trans.net$itimes[-1]))))
transmission.events.df <- data.frame(number.trans = annual.infections,
                                     calendaryear = year.bins)

# Timing of internal nodes
aged.dater.tree.125.random <- picante::node.age(dater.tree.125.random)
annual.internal.nodes.125.random <- as.numeric(table(floor(aged.dater.tree.125.random$Ti)))
year.bins.intnodes.125.random <- as.numeric(names(table(floor(aged.dater.tree.125.random$Ti))))
internal.nodes.125.random.df <- data.frame(number.nodes = annual.internal.nodes.125.random,
                                          calendaryear = year.bins.intnodes.125.random)

years.df <- data.frame(calendaryear = 1985:2017)
trans.df <- dplyr::left_join(x = years.df,
                             y = transmission.events.df) %>%
  replace_na(list(number.trans = 0))
trans.and.nodes.125.random.df <- dplyr::left_join(x = trans.df,
                                                 y = internal.nodes.125.random.df) %>%
  replace_na(list(number.nodes = 0))


trans.and.nodes.long.125.random.df <- gather(trans.and.nodes.125.random.df,
                                            key = "Events",
                                            value = "Number",
                                            number.trans:number.nodes,
                                            factor_key = TRUE)


transandnodes.125.random.plot <- ggplot(data = trans.and.nodes.long.125.random.df,
                                       aes(x = calendaryear,
                                           y = Number,
                                           colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time")
print(transandnodes.125.random.plot)

ggsave(filename = "transandnodes.125.random.pdf",
       plot = transandnodes.125.random.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")

## Adding density to dataset trans.and.nodes.long.df
trans.and.nodes.long.125.random.enriched.df <- trans.and.nodes.long.125.random.df %>%
  group_by(Events) %>%
  mutate(total.events = sum(Number),
         percentage.events = Number / total.events)

# New plot comparing densities  
transandnodes.fraction.125.random.plot <- ggplot(data = dplyr::filter(trans.and.nodes.long.125.random.enriched.df,
                                                                     calendaryear >= 1985,
                                                                     calendaryear <= 2017),
                                                aes(x = calendaryear,
                                                    y = percentage.events,
                                                    colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Transmission events",
                                "Internal nodes")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2017),
                     breaks = seq(from = 1985,
                                  to = 2017,
                                  by = 5)) +
  xlab("Time") +
  ylab("Fraction")
print(transandnodes.fraction.125.random.plot)

ggsave(filename = "transandnodesfraction125random.pdf",
       plot = transandnodes.fraction.125.random.plot,
       path = paste0(getwd(), "/plots"),
       width = 12, height = 12, units = "cm")




# RMSE of percentages
RMSE.full <- sqrt(sum((trans.and.nodes.long.enriched.df$percentage.events[1:31] - trans.and.nodes.long.enriched.df$percentage.events[32:62])^2) / 31)

# RMSE of percentages
RMSE.50 <- sqrt(sum((trans.and.nodes.long.50.random.enriched.df$percentage.events[1:31] - trans.and.nodes.long.50.random.enriched.df$percentage.events[32:62])^2) / 31)

# RMSE of percentages
RMSE.25 <- sqrt(sum((trans.and.nodes.long.25.random.enriched.df$percentage.events[1:31] - trans.and.nodes.long.25.random.enriched.df$percentage.events[32:62])^2) / 31)

# RMSE of percentages
RMSE.125 <- sqrt(sum((trans.and.nodes.long.125.random.enriched.df$percentage.events[1:31] - trans.and.nodes.long.125.random.enriched.df$percentage.events[32:62])^2) / 31)


####################################
## Updated figure 5C
####################################
trans.and.nodes.long.enriched.df$Scenario <- "Transmission events"
trans.and.nodes.long.enriched.df$Scenario[trans.and.nodes.long.enriched.df$Events == "number.nodes"] <- "100% sequence coverage"
trans.and.nodes.long.50.random.enriched.df$Scenario <- "Random 50% sequence coverage"
trans.and.nodes.long.25.random.enriched.df$Scenario <- "Random 25% sequence coverage"
trans.and.nodes.long.125.random.enriched.df$Scenario <- "Random 12.5% sequence coverage"

trans.and.nodes.long.enriched.combined <- rbind(trans.and.nodes.long.enriched.df,
                                                dplyr::filter(trans.and.nodes.long.50.random.enriched.df,
                                                              Events == "number.nodes"),
                                                dplyr::filter(trans.and.nodes.long.25.random.enriched.df,
                                                              Events == "number.nodes"),
                                                dplyr::filter(trans.and.nodes.long.125.random.enriched.df,
                                                              Events == "number.nodes"))


trans.and.nodes.long.enriched.combined$Scenario <- factor(trans.and.nodes.long.enriched.combined$Scenario,
                                                          levels = rev(c("Transmission events",
                                                                     "100% sequence coverage",
                                                                     "Random 50% sequence coverage",
                                                                     "Random 25% sequence coverage",
                                                                     "Random 12.5% sequence coverage")))


#                     Red      Dark blue   Blue      Light blue  Baby blue
custom.colours <- c("#C90001", "#00427A", "#096CBF", "#4FAFFF", "#A3D5FF")

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
        legend.title = element_blank()) +
  scale_color_manual(values = c("Transmission events" = custom.colours[1],
                                "100% sequence coverage" = custom.colours[2],
                                "Random 50% sequence coverage" = custom.colours[3],
                                "Random 25% sequence coverage" = custom.colours[4],
                                "Random 12.5% sequence coverage" = custom.colours[5])) 
print(transandnodescombined.plot)

ggsave(filename = "transandnodescombined.pdf",
       plot = transandnodescombined.plot,
       path = paste0(getwd(), "/plots"),
       width = 20, height = 12, units = "cm")
