require(ape)
require(vegan)
require(picante)
require(ade4)
require(bipartite)
require(ggplot2)
require(reshape2)
require(Deducer)

# Format Network Matrix

Associations <- read.table("All_Associations.txt",header=TRUE,sep="\t")
IGV <- read.table("Potential_Intragenomic_Variants.txt",sep="\t")
    Associations_Keep <- droplevels(Associations[! Associations$Phylotype  %in% IGV[,1],])
    Associations_Drop <- droplevels(Associations[Associations$Phylotype  %in% IGV[,1],])


    Coral_Matrix <- as.data.frame.matrix(table(Associations_Drop[,2:3]))



    Coral_Matrix <- Coral_Matrix[rowSums(Coral_Matrix) > 9,]
      Coral_Matrix <- Coral_Matrix[,colSums(Coral_Matrix)  > 2] # At least 3 records
      Coral_Matrix_Percent <- Coral_Matrix/rowSums(Coral_Matrix)

# Calculate Modules

  nreps <- 1000
    AllMods <- list()


NumMods <- function(Mods) {
  
  ModuleData <- Mods@modules[2:nrow(Mods@modules),3:ncol(Mods@modules)]
  return(nrow(ModuleData))
  
}  

for (i in 1:nreps) {
  print(i)
  AllMods[[i]] <- computeModules(Coral_Matrix_Percent)
  print(NumMods(AllMods[[i]]))
  
}

Mods <- AllMods[[1]]


Output <- matrix(0,nrow=nreps,ncol=3)

for(i in 1:nreps) {
  
  print(i)
  
  Output[i,1] <- paste("Rep",i,sep="")
  
  Output[i,2] <- as.numeric(as.character(NumMods(AllMods[[i]])))
  
  Output[i,3] <- as.numeric(as.character(AllMods[[i]]@likelihood))
  
}

Output <- as.data.frame(Output)
colnames(Output) <- c("Rep","NumMods","Q")
Output$NumMods <- as.numeric(as.character(Output$NumMods))
Output$Q <- as.numeric(as.character(Output$Q))



ggplot(Output,aes(x=NumMods,y=Q)) + geom_point() + theme_bw()




######################


Mods <- AllMods[[1]]


    plotModuleWeb(Mods)
    Mods@likelihood

    ModuleData <- Mods@modules[2:nrow(Mods@modules),3:ncol(Mods@modules)]
      row.names(ModuleData) <- paste("Mod",seq(1:nrow(ModuleData)),sep="")
      colnames(ModuleData) <- c(row.names(Coral_Matrix_Percent),colnames(Coral_Matrix_Percent))

    CoralMods     <- ModuleData[,1:nrow(Coral_Matrix_Percent)]
    PhylotypeMods <- ModuleData[,(1 + nrow(Coral_Matrix_Percent)):ncol(ModuleData)]

    ModuleSizes <- rowSums(ModuleData != 0)
    ModuleSizes_Coral <- rowSums(CoralMods != 0)
    ModuleSizes_Phylotype <- rowSums(PhylotypeMods != 0)

Output <- matrix(0,nrow=ncol(ModuleData),ncol=3)

for(i in 1:nrow(Output)) {
  
  Output[i,1] <- colnames(ModuleData)[i]
  
  if(colnames(ModuleData)[i] %in% colnames(CoralMods)) {Output[i,2] <- "Coral"}
  else if(colnames(ModuleData)[i] %in% colnames(PhylotypeMods)) {Output[i,2] <- "Symbiodinium"}
  
  NodeData <- ModuleData[,i]
  
  Output[i,3] <- as.character(names(NodeData[NodeData != 0]))
  
}

ModuleAssignments <- as.data.frame(Output)
colnames(ModuleAssignments) <- c("Node","Type","Module")

Metadata <- read.table("NetworkMetadata.txt",header=TRUE,sep="\t")

ModuleAssignments$BRI     <- Metadata$Coral_BRI[match(ModuleAssignments$Node,Metadata$Node)]
ModuleAssignments$TT      <- Metadata$Phylotype_TT[match(ModuleAssignments$Node,Metadata$Node)]
ModuleAssignments$Records <- Metadata$Coral_Number_Records[match(ModuleAssignments$Node,Metadata$Node)]

ModulesKeep <- names(ModuleSizes[ModuleSizes > 3])
ModuleAssignmentsKeep <- ModuleAssignments[ModuleAssignments$Module %in% ModulesKeep,]

### PLOT MODULE TT AGAINST MODULE BRI ###

AllMods <- unique(ModuleAssignments$Module)

Output <- matrix(0,nrow=length(unique(ModuleAssignments$Module)),ncol=9)
colnames(Output) <- c("Module","Num_Corals","Num_Phylotypes","Num_Phylotypes_with_TT","Num_Phylotypes_without_TT","Percent_Observations_with_TT","Weighted_Mean_TT","Mean_BRI","Weighted_Mean_BRI")

for (i in 1:length(AllMods)) {
  
  Output[i,1] <- as.character(AllMods[i]) # Module number
  
  Mod_i_Members <- ModuleAssignments[ModuleAssignments$Module == AllMods[i],] 
  
  Coral_Matrix_Mod <- Coral_Matrix[row.names(Coral_Matrix) %in% Mod_i_Members$Node,colnames(Coral_Matrix) %in% Mod_i_Members$Node,drop=FALSE]
  
  Output[i,2] <- as.numeric(as.character(nrow(Coral_Matrix_Mod))) # Number of corals in mod
  
  Output[i,3] <- as.numeric(as.character(ncol(Coral_Matrix_Mod))) # Number of phylotypes in mod
  
  Mod_PhylotypeCounts <- as.data.frame(colSums(Coral_Matrix_Mod))
  colnames(Mod_PhylotypeCounts) <- "Count"
  Mod_PhylotypeCounts$TT <- Metadata$Phylotype_TT[match(row.names(Mod_PhylotypeCounts),Metadata$Node)]
  
  TT_Negative <- Mod_PhylotypeCounts[is.na(Mod_PhylotypeCounts$TT),]
  TT_Positive <- Mod_PhylotypeCounts[complete.cases(Mod_PhylotypeCounts)==TRUE,]
  
  Output[i,4] <- nrow(TT_Positive) # Number of phylotypes with TT
  Output[i,5] <- nrow(TT_Negative) # Number of phylotypes without TT
  
  PercentWithTT <- sum(TT_Positive$Count)/sum(Mod_PhylotypeCounts$Count)
  
  Output[i,6] <- PercentWithTT # Percent of observation in module with phylotypes with known TT
  
  if (PercentWithTT == 0) {Output[i,7] <- NA}
  else if (PercentWithTT > 0) {
    TT_Positive$Percent <- TT_Positive$Count/sum(TT_Positive$Count)
    Weighted_Mean_TT <- sum(TT_Positive$TT*TT_Positive$Percent)
    
    Output[i,7] <- Weighted_Mean_TT # Weighted mean TT in module
  }
  
  Mod_CoralCounts <- as.data.frame(rowSums(Coral_Matrix_Mod))
  colnames(Mod_CoralCounts) <- "Count"
  
  Mod_CoralCounts$BRI <- Metadata$Coral_BRI[match(row.names(Mod_CoralCounts),Metadata$Node)]
  
  Output[i,8] <-  mean(Mod_CoralCounts$BRI) # Mean BRI
  
  Mod_CoralCounts$Percent <- Mod_CoralCounts$Count/sum(Mod_CoralCounts$Count)
  
  Output[i,9] <- sum(Mod_CoralCounts$BRI*Mod_CoralCounts$Percent) # Weighted Mean BRI
  
  
}

Output <- as.data.frame(Output)

ModData <- Output[complete.cases(Output),]
for (i in 2:ncol(ModData)) {
  ModData[,i] <- as.numeric(as.character(ModData[,i]))
}

ggplot(ModData,aes(x=Weighted_Mean_TT,y=Mean_BRI)) + geom_point() + geom_smooth(method="lm")

### PHYLO MODULE CORRELATION ###

SymTree <- read.tree("Tree_AllPhylotypes.tre")
SymTreeDistances <- cophenetic.phylo(SymTree)
SymTreeDistancesMelt <- melt(SymTreeDistances)
colnames(SymTreeDistancesMelt) <- c("Phy1","Phy2","Distance")

SymDistance <- SymTreeDistancesMelt[SymTreeDistancesMelt$Phy1 %in% ModuleAssignments$Node & SymTreeDistancesMelt$Phy2 %in% ModuleAssignments$Node,]

SymDistance$Phy1_Mod <- ModuleAssignments$Module[match(SymDistance$Phy1,ModuleAssignments$Node)] 
SymDistance$Phy2_Mod <- ModuleAssignments$Module[match(SymDistance$Phy2,ModuleAssignments$Node)] 
SymDistance$Comp <- 0

for(i in 1:nrow(SymDistance))
{
  if(SymDistance$Phy1_Mod[i] == SymDistance$Phy2_Mod[i]) {SymDistance$Comp[i] <- "Within"}
  else if(SymDistance$Phy1_Mod[i] != SymDistance$Phy2_Mod[i]) {SymDistance$Comp[i] <- "Between"}
}

CoralTree <- read.nexus("CoralShortPhylo320.nex")  
CoralTreeDistances <- cophenetic.phylo(CoralTree)
CoralTreeDistancesMelt <- melt(CoralTreeDistances)
colnames(CoralTreeDistancesMelt) <- c("Coral1","Coral2","Distance")
CoralTreeDistancesMelt$Coral1 <- gsub("_"," ",CoralTreeDistancesMelt$Coral1)
CoralTreeDistancesMelt$Coral2 <- gsub("_"," ",CoralTreeDistancesMelt$Coral2)

CoralDistance <- CoralTreeDistancesMelt[CoralTreeDistancesMelt$Coral1 %in% ModuleAssignments$Node & CoralTreeDistancesMelt$Coral2 %in% ModuleAssignments$Node,]

CoralDistance$Coral1_Mod <- ModuleAssignments$Module[match(CoralDistance$Coral1,ModuleAssignments$Node)] 
CoralDistance$Coral2_Mod <- ModuleAssignments$Module[match(CoralDistance$Coral2,ModuleAssignments$Node)] 
CoralDistance$Comp <- 0

for(i in 1:nrow(CoralDistance))
{
  if(CoralDistance$Coral1_Mod[i] == CoralDistance$Coral2_Mod[i]) {CoralDistance$Comp[i] <- "Within"}
  else if(CoralDistance$Coral1_Mod[i] != CoralDistance$Coral2_Mod[i]) {CoralDistance$Comp[i] <- "Between"} 
}

CoralDistance <- CoralDistance[CoralDistance$Distance != 0,]
SymDistance <- SymDistance[SymDistance$Distance != 0,]

P1 <- ggplot(CoralDistance,aes(x=Distance)) + geom_density(aes(col=Comp,fill=Comp),alpha=0.3) + theme_bw() + scale_fill_manual(values=c("#67a9cf","#ef8a62")) + scale_color_manual(values=c("#67a9cf","#ef8a62")) + labs(x="Phylogenetic Distance between Coral Species")
P2 <- ggplot(SymDistance,aes(x=Distance)) + geom_density(aes(col=Comp,fill=Comp),alpha=0.3) + theme_bw() + scale_fill_manual(values=c("#67a9cf","#ef8a62")) + scale_color_manual(values=c("#67a9cf","#ef8a62")) + labs(x="Phylogenetic Distance between Symbiodinium Phylotypes")

Coral_Within <- CoralDistance[CoralDistance$Comp == "Within",]
Coral_Between <- CoralDistance[CoralDistance$Comp == "Between",]

Sym_Within <- SymDistance[SymDistance$Comp == "Within",]
Sym_Between <- SymDistance[SymDistance$Comp == "Between",]

# Perm T-Tests

CoralTTest <- t.test(Coral_Within$Distance,Coral_Between$Distance)
SymTTest   <- t.test(Sym_Within$Distance,Sym_Between$Distance)

#Symbiodinium
nreps <- 10000
TTest_Output <- matrix(0,nrow=nreps,ncol=3)
TTest_Output <- as.data.frame(TTest_Output)
colnames(TTest_Output) <- c("Run","TestStat","p")
Rows <- 1:18360

for (i in 1:nreps) {
  
  if(i%%1000 == 0) {print(i)}
  
  WithinRows <- sample.int(18360, 2858)
  BetweenRows <- setdiff(Rows,WithinRows)
  
  WithinData <- SymDistance[WithinRows,]
  BetweenData <- SymDistance[BetweenRows,]
  
  ttest <- t.test(WithinData$Distance,BetweenData$Distance)
  
  TTest_Output[i,1] <- as.character(i)
  TTest_Output[i,2] <- ttest$statistic
  TTest_Output[i,3] <- ttest$p.value
  
}

TTest_Output$TestStat <- as.numeric(as.character(TTest_Output$TestStat))
TTest_Output$p <- as.numeric(as.character(TTest_Output$p))s

P3 <- ggplot(TTest_Output,aes(x=TestStat)) + geom_density(fill="#999999",alpha=0.3) + geom_point(aes(x=SymTTest$statistic,y=0),col="#b2182b",pch=10,size=6) + theme_bw() + labs(x="Test Statistic")

#Coral

nreps <- 10000
TTest_Output_C <- matrix(0,nrow=nreps,ncol=3)
TTest_Output_C <- as.data.frame(TTest_Output_C)
colnames(TTest_Output_C) <- c("Run","TestStat","p")
Rows <- 1:21462

for (i in 1:nreps) {
  
  if(i%%1000 == 0) {print(i)}
  
  WithinRows <- sample.int(21462, 3648)
  BetweenRows <- setdiff(Rows,WithinRows)
  
  WithinData <- SymDistance[WithinRows,]
  BetweenData <- SymDistance[BetweenRows,]
  
  ttest <- t.test(WithinData$Distance,BetweenData$Distance)
  
  TTest_Output_C[i,1] <- as.character(i)
  TTest_Output_C[i,2] <- ttest$statistic
  TTest_Output_C[i,3] <- ttest$p.value
  
}

P4 <- ggplot(TTest_Output_C,aes(x=TestStat)) + geom_density(fill="#999999",alpha=0.3) + geom_point(aes(x=CoralTTest$statistic,y=0),col="#b2182b",pch=10,size=6) + theme_bw() + labs(x="Test Statistic")

grid.arrange(P1,P4,P2,P3,nrow=2)

################
TT_Mean <- aggregate(TT ~ Module, ModuleAssignmentsKeep, function(x) mean = mean(x))

ModuleAssignmentsKeep2 <- ModuleAssignmentsKeep[ModuleAssignmentsKeep$Module %in% TT_Mean$Module,] 

TT_SEM  <- aggregate(TT ~ Module, ModuleAssignmentsKeep2, function(x) se = sd(x)/sqrt(length(x)))
BRI_Mean <- aggregate(BRI ~ Module, ModuleAssignmentsKeep2, function(x) mean = mean(x))
BRI_SEM  <- aggregate(BRI ~ Module, ModuleAssignmentsKeep2, function(x) se = sd(x)/sqrt(length(x)))
Length <- aggregate(Node ~ Module, ModuleAssignmentsKeep2, function(x) length = (length(x)))

Merged <- cbind(BRI_Mean,BRI_SEM,TT_Mean,TT_SEM,Length)[,c(1,2,4,6,8,10)]
colnames(Merged) <- c("Module","BRI_Mean","BRI_SEM","TT_Mean","TT_SEM","Size")
Merged[,2] <- as.numeric(as.character(Merged[,2]))
Merged[,3] <- as.numeric(as.character(Merged[,3]))
Merged[,4] <- as.numeric(as.character(Merged[,4]))
Merged[,5] <- as.numeric(as.character(Merged[,5]))

ggplot(Merged,aes(x=BRI_Mean,y=TT_Mean)) + geom_point(aes(size=Size),col="#252525") + geom_errorbar(aes(x=BRI_Mean,ymin=TT_Mean-TT_SEM,ymax=TT_Mean+TT_SEM),alpha=0.7,col="#252525") + geom_smooth(method="lm",col="orange",fill="orange",alpha=0.1) +  geom_errorbarh(aes(y=TT_Mean,xmin=BRI_Mean-BRI_SEM,xmax=BRI_Mean+BRI_SEM),alpha=0.7) + labs(x="Mean Module BRI",y="Mean Module TT")  + theme_bw() + theme(legend.position="none") + geom_text(aes(label=Module),col="orangered",fontface=2)

ggplot(ModuleAssignments,aes(x=Module,y=BRI)) + geom_boxplot()


