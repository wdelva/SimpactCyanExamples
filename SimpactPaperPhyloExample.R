# Before running the code below, make sure you have followed the instructions for installing all the required software,
# as explained in the README file.


setwd("/path/to/your/working_directory/") # Change this to your own working directory

# Make sure compiled tools (Seq-Gen and FastTree) are in same working directory

## Load required packages
library(RSimpactCyan)
library(RSimpactHelper)
library(phangorn)
library(treedater)
library(picante)
library(igraph)
library(geomnet)
library(lubridate)
library(ggtree)
library(ggplot2)
library(tidyr)


inputvector <- c(123, -0.52, -0.05, 5, 7, 3, 0.25, -0.3, -0.1, 
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) 



#######################
# Step 1: Run Simpact #
#######################

## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
#
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 50, 
                                 population.nummen = 2000, 
                                 population.numwomen = 2000,
                                 hivseed.time = 10, 
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, 
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 debut.debutage = 15
)

# # Assumption of nature of sexual network
# #########################################
#
cfg.list["population.msm"] = "no"


# # Sexual behaviour
# ###################
#
seedid <- inputvector[1]

cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)


# # HIV transmission
# ###################
#

cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)

# Disease progression > may be remove in parameter to estimates

cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)


#
# # Demographic
# ##############
#

cfg.list["conception.alpha_base"] <- inputvector[18]


# # Assumptions to avoid negative branch lengths
# ###############################################
# # + sampling == start ART
# # when someone start ART, he/she is sampled and becomes non-infectious

cfg.list["monitoring.fraction.log_viralload"] <- 0


#
# ## Add-ons
#
### BEGIN Add-on
cfg.list["formation.hazard.agegapry.baseline"] <- 2
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
cfg.list["dropout.interval.dist.type"] <- "uniform"
cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
cfg.list["diagnosis.baseline"] <- -2


cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23
cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23
cfg.list["person.eagerness.man.dist.gamma.b"] <- 45
cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45

#### END Add-ons


# # ART intervention
# ###################
#
# # ART acceptability parameter and the ART interventions

cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6

# Let's introduce ART
art.intro <- list()
art.intro["time"] <- 20
art.intro["diagnosis.baseline"] <- -2
art.intro["monitoring.cd4.threshold"] <- 100

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis
# Gradual increase in CD4 threshold. in 2005:200. in 2010:350. in 2013:500
art.intro1 <- list()
art.intro1["time"] <- 22
art.intro1["diagnosis.baseline"] <- -2
art.intro1["monitoring.cd4.threshold"] <- 150

art.intro2 <- list()
art.intro2["time"] <- 25
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 30
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 33
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 36
art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access

intervention <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)

# Events
cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3


#
# # # # Run Simpact
results <- simpact.run(configParams = cfg.list,
                       destDir = "temp",
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)

datalist <- readthedata(results)

###########################################
# Step 2: Construct transmission networks #
###########################################


# Produce a list of transmission networks in epi object format

simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)


# Check which transmission network with at least 3 individuals
# and return the transmission object

smallest.branches <- rep(NA, times = length(simpact.trans.net))
for (list.element in 1:length(simpact.trans.net)){
  net.list <- simpact.trans.net[[list.element]]
  if(length(net.list$id) > 2){
    tree.tr <- epi2tree2(net.list)
    smallest.branch <- min(tree.tr$edge.length)
    smallest.branches[list.element] <- smallest.branch
  }
}




###############################
# Step 3: Sequence simulation #
###############################


# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed <-  123

trans.net <- simpact.trans.net # all transmission networks

num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network

# constrained to rename IDs to -1, 0, 1, 2, ...

num.i <- vector() # i_th seed in the list of seeds

for(i in 1:length(trans.net)){
  
  tree.n <- trans.net[[i]] # transmission network for i^th seed
  
  if(nrow(as.data.frame(tree.n)) >= 50){ # consider transmission networks with at least 50 individuals
    tree.i <- epi2tree2(tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
    num.i <- c(num.i,i)
    
    # Save the transmission tree
    write.tree(tree.i, file = paste("tree.model1.seed",i,".nwk", sep = ""))
    
    tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    
    seq.rand <- 1
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)
    
    n.tr <- 1
    
    seed.id <- tree.n$id[1]
    
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))
    
    system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > B.EpidemicSequences_seed_",i,".fasta",sep = ""))
    
    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # s: scale which is the substitution rate of pol gene
    # z: seed
    
    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network
    
    write.csv(id.samplingtime,file=paste("samplingtimes_seed_",i,".csv", sep = ""))
    
  }
}

# Chosen transmission networks with at least 50 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds

# > IDs.transm
# [1]  11 15
# > nrow(data.frame(simpact.trans.net[[11]]))
# [1] 426
# > nrow(data.frame(simpact.trans.net[[15]]))
# [1] 52



#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################


# 4.1. Build phylogenetic trees

for (j in 1:length(IDs.transm)){
  
  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log
  
  # system(paste("./iqtree-omp -s", paste("B.EpidemicSequences_seed_",id.trans,".fasta", sep = ""), " -nt AUTO -alrt 1000 -bb 1000"))
  
  # Compiling FastTree
  # gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
  
  system(paste("./FastTree -gtr -nt <", paste("B.EpidemicSequences_seed_",id.trans,".fasta", sep = ""), paste(">B.EpidemicSequences_seed_",id.trans,".fasta.tree", sep = "")))
  
  
}


# 4.2. Calibrated internal nodes of the phylogenetic trees

# Function to tranform sampling time into calendar time - dates in named vector to be handled
# by treedater during internal nodes calibration

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1977+40-dates$V2 # dates datalist$itable$population.simtime[1] - dates$V2 + 1977
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}



# For ALL transmission networks IDs.transm #
############################################


# Not necessary for the purpose of this example to prove the concept
# of relatedness of HIV phylodynamics and transmission network

for (j in 1:length(IDs.transm)){
  
  id.trans <- IDs.transm[j]
  
  tree.const <- read.tree(paste("B.EpidemicSequences_seed_",id.trans,".fasta.tree", sep = ""))
  
  samp.dates <- read.csv(paste("samplingtimes_seed_",id.trans,".csv", sep = ""))
  
  time.samp <- dates.Transform.NamedVector(dates=samp.dates)
  
  tree.tips <- as.numeric(tree.const$tip.label)
  
  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }
  
  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012) # s is the length of sequence
  
  save(dater.tree, file = paste("dated.tree.object_seed_",id.trans,".Rdata", sep = ""))
  
  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  
}


# For ONE transmission network IDs.transm[1] #
##############################################

# For the purpose of this exercise we will consider ONE transmission network 

# Calibration of the phylogenetic tree produced with sequences from
# transmission network 11 (IDs.transm[1])

id.trans <- IDs.transm[1] # 11

tree.const <- read.tree(paste("B.EpidemicSequences_seed_",id.trans,".fasta.tree", sep = ""))

samp.dates <- read.csv(paste("samplingtimes_seed_",id.trans,".csv", sep = ""))

time.samp <- dates.Transform.NamedVector(dates=samp.dates)

tree.tips <- as.numeric(tree.const$tip.label)


# Sort sampling dates obtained with dates.Transform.NamedVector function
# in same order as in the phylogenetic tree

Ord.tree.dates <- vector()

for(i in 1:length(tree.tips)){
  for(j in 1:length(time.samp)){
    if(tree.tips[i] == samp.dates$V1[j]){
      Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
    }
  }
}




# Use of library(treedater) to calibrate internal nodes

dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012) # s is the length of sequence

# dater.tree is an object with much information on the tree and evolutionary parameters

# Save the object
save(dater.tree, file = "dater.tree.RData")

# We can save the tree
write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))

# Node age with picante package

N <- node.age(dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date



#####################################################################################
# Step 5: Compute transmission events and internal nodes in one year time intervals #
#####################################################################################


# THis help to assess transmission events and internal nodes distribution in one year time intervals

# We choose  transmission network 11 which is IDs.transm[1]


## Transmission network of seed 11
trans.net <- simpact.trans.net

tra.net.11 <- trans.net[[11]]

tra.net.11$dtimes <- abs(tra.net.11$dtimes-40)+1977
tra.net.11$itimes <- abs(tra.net.11$itimes-40)+1977

min.val = 1977
max.val = round(max(tra.net.11$itimes))



# (i) NUmber of internal nodes & transmission events

step.int=1
d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(tra.net.11)
dt.node.age.dt <- int.node.age

numb.tra <- vector() # initialize transmission events
i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  inf <- 1976+i
  sup <- 1977+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes <= sup & dat.f.trans$itimes  > inf),]
  numb.i <- nrow(dat.f.trans.i)
  numb.tra <- c(numb.tra, numb.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}


# (ii) Transmission network

graph.build <- as.data.frame(trans.net[[11]])

graph.build[,4] <- as.character(graph.build$parent) # donors
graph.build[,3] <- as.character(graph.build$id) # recipients
gag = as.matrix(graph.build)
gag = gag[-1,] # remove universall seed -1
ga.graph = graph.edgelist(gag[,4:3])

V(ga.graph)$color <- "red"

transNet.yrs.Old <- ga.graph




#########################################################################
# Step 6: Creating the panels A, B and C of the figure for this example #
#########################################################################


# Objects needed for the figure:
SimpactPaperPhyloExample <- list()
SimpactPaperPhyloExample$transNet.yrs.Ord <- transNet.yrs.Old
SimpactPaperPhyloExample$dater.tree <- dater.tree
SimpactPaperPhyloExample$i.vec <- i.vec
SimpactPaperPhyloExample$int.node.vec <- int.node.vec
SimpactPaperPhyloExample$numb.tra <- numb.tra
save(SimpactPaperPhyloExample, file = "SimpactPaperPhyloExample.RData")

load(file = "/path/to/your/working_directory/SimpactPaperPhyloExample.RData")
# A. Transmission network

network <- SimpactPaperPhyloExample$transNet.yrs.Ord
edges <- as.data.frame(get.edgelist(network))
names(edges) <- c("infector", "infectee")
vertices <- as.data.frame(vertex_attr(network))
network.df <- fortify(as.edgedf(edges), vertices)

transmissionnetwork.plot <- ggplot(data = network.df,
                                   aes(from_id = from_id,
                                       to_id = to_id)) +
  geom_net(directed = TRUE,
           size = 2.5,
           layout.alg = "kamadakawai",
           layout.par = list(niter = 2000,
                             sigma = 200,
                             kkconst = 10),
           alpha = 1,
           arrowsize = 0.7,
           arrowgap = 0.005,
           ecolour = "darkgrey",
           colour = "black") +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  ylab("")
print(transmissionnetwork.plot)

ggsave(filename = "network_vsc.pdf",
       plot = transmissionnetwork.plot,
       path = "/path/to/your/working_directory/plots",
       width = 20, height = 30, units = "cm")



# B. Phylogenetic tree

tree <- SimpactPaperPhyloExample$dater.tree
class(tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
sim.start.year <- 1987
first.transmission <- min(SimpactPaperPhyloExample$dater.tree$sts)
mrsd <- max(SimpactPaperPhyloExample$dater.tree$sts)

dates <- format(date_decimal(c(mrsd, first.transmission)), "%Y-%m-%d")
tree$root.edge <-  - sim.start.year
phylotree.plot <- ggtree(tree,
                         mrsd = dates[1],
                         size = 0.05) + 
  theme_tree2() +
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

ggsave(filename = "tree_vsc.pdf",
       plot = phylotree.plot,
       path = "/path/to/your/working_directory/plots",
       width = 10, height = 15, units = "cm")




# C. Transmission event versus internal nodes

calendaryear.isrelevant <- SimpactPaperPhyloExample$i.vec >= 1987
calendaryear <- SimpactPaperPhyloExample$i.vec[calendaryear.isrelevant]
intern.nodes <- SimpactPaperPhyloExample$int.node.vec[calendaryear.isrelevant]
trans.events <- SimpactPaperPhyloExample$numb.tra[calendaryear.isrelevant]

trans.and.nodes.df <- data.frame(calendaryear = calendaryear,
                                 intern.nodes = intern.nodes,
                                 trans.events = trans.events)
trans.and.nodes.long.df <- gather(trans.and.nodes.df,
                                  key = "Events",
                                  value = "Number",
                                  intern.nodes:trans.events,
                                  factor_key = TRUE)


transandnodes.plot <- ggplot(data = trans.and.nodes.long.df,
                             aes(x = calendaryear,
                                 y = Number,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Internal nodes",
                                "Transmission events")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) +
  xlab("Time")
print(transandnodes.plot)

ggsave(filename = "events_vsc.pdf",
       plot = transandnodes.plot,
       path = "/path/to/your/working_directory/plots",
       width = 10, height = 15, units = "cm")
