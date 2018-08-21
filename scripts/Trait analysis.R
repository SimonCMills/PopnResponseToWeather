## PGLS analysis, relating contribution of weather scores to traits and phylogeny
##
##  Phylogeny from "Roquet et al. (2014) One tree to rule them all: a phylogenetic
##  dataset for the European Tetrapoda"

# Housekeeping
library(dplyr); library(caper); library(xlsx)

## Data ----
# get species key from full dataset
speciesKey <- readRDS("files/birdDF.rds") %>%
    ungroup %>%
    dplyr::select(species, species_LB) %>%
    unique
# read in tree
tree <- read.nexus("files/phylo_treebase.nex")
# read in formatted trait database
traitDB <- readRDS("files/traitDB.rds")

# read in model outputs, for full variable set and reduced, respectively
dataset_full <- readRDS("files/model outputs/fullModels_explanatoryCapacity.rds") %>%
    left_join(., speciesKey) %>%
    inner_join(., traitDB)

dataset_red <-  readRDS("files/model outputs/redModels_explanatoryCapacity.rds") %>%
    left_join(., speciesKey) %>%
    inner_join(., traitDB)

## format tree and model outputs for pgls ----
tree_sp <- tree$tip.label
# species in dataset
data_sp <- dataset_red$species_LB %>%
    gsub(" ", "_", .)

## All species bar hooded crow and lesser redpoll are in the phylogeny, but 10 
## have different LB which needs correcting and formatting
missingPhylo <- dataset_red$species[!data_sp %in% tree_sp]
(updatedID <- data_frame(species = missingPhylo,
                        LB = data_sp[!data_sp %in% tree_sp],
                        phylo = c("Miliaria_calandra", "Parus_montanus", "Larus_graellsii",
                                  "Parus_cristatus", "Saxicola_torquata", "Parus_palustris",
                                  "Larus_ridibundus", "Delichon_urbica", "Parus_ater", 
                                  "Parus_caeruleus")))


# format full variable set df
dataset_full_pgls <- left_join(dataset_full, updatedID) %>%
    dplyr::select(-LB) %>%
    mutate(phylo = ifelse(is.na(phylo), gsub(" ", "_", species_LB), phylo)) %>%
    as.data.frame %>%
    filter(complete.cases(.))
row.names(dataset_full_pgls) <- dataset_full_pgls$phylo

# check: no species missing info in final df (TRUE)
dataset_full_pgls %>% filter(!complete.cases(.))

# format reduced variable set df
dataset_red_pgls <- left_join(dataset_red, updatedID) %>%
    dplyr::select(-LB) %>%
    mutate(phylo = ifelse(is.na(phylo), gsub(" ", "_", species_LB), phylo)) %>%
    as.data.frame
row.names(dataset_red_pgls) <- dataset_red_pgls$phylo

# check: no species missing info in final df (TRUE)
dataset_red_pgls %>% filter(!complete.cases(.))

## check: all species in dataset now have matching phylogeny
data_updated_sp <- dataset_red_pgls$phylo 
data_updated_sp[!data_updated_sp %in% tree_sp]

## prune tree to species in respective datasets
toDrop_full <- which(!tree$tip.label %in% dataset_full_pgls$phylo )
tree_full_pruned <- drop.tip(tree, toDrop_full)

toDrop_red <- which(!tree$tip.label %in% dataset_red_pgls$phylo )
tree_red_pruned <- drop.tip(tree, toDrop_red)

## combine to comparative data object
redVars_comparative <- comparative.data(tree_red_pruned, dataset_red_pgls, names.col = phylo)
fullVars_comparative <- comparative.data(tree_full_pruned, dataset_full_pgls, names.col = phylo)

## run PGLS models ----
# formula
form <- 1-var_quotient ~ log(body_mass) + multibrood + LD_migrant + life_span +
    diet_vert + diet_invert + diet_plant + hab_forest + hab_open + hab_aquatic

f_redVars <- pgls(form, redVars_comparative, lambda="ML")
f_fullVars <- pgls(form, fullVars_comparative, lambda="ML")

summary(f_redVars)
summary(f_fullVars)

prof <- pgls.profile(f_fullVars, "lambda")
plot(prof)

prof2 <- pgls.profile(f_redVars, "lambda")
plot(prof2)

par(mfrow=c(2,2)); plot(f_fullVars)
par(mfrow=c(2,2)); plot(f_redVars)

# tabulate and save to excel
tabulateSummary <- function(model) {
    summ <- summary(model)
    lambda <- summ$param.CI$lambda
    lambda_df <- data_frame(Coefficient = "lambda", 
               Estimate = round(lambda$opt, 3),
               uprLwr = paste0(round(lambda$ci.val, 3), collapse=", "), 
               CI = paste0("(", uprLwr, ")"), 
               pVal = NA) %>%
        dplyr::select(-uprLwr)
    R2 <- data_frame(Coefficient = "R2", Estimate=summ$r.squared, CI = NA, pVal=NA)
    coef(summ) %>%
        as.data.frame %>%
        mutate(Coefficient = row.names(.), 
               lwr = round(Estimate - 1.96*`Std. Error`, 3), 
               upr = round(Estimate + 1.96*`Std. Error`, 3),
               Estimate = round(Estimate, 3), 
               CI = paste0("(", lwr, ", ", upr, ")"),
               pVal = round(`Pr(>|t|)`, 3)) %>%
        dplyr::select(Coefficient, Estimate, CI, pVal) %>%
        bind_rows(., lambda_df, R2)
} 

# write table ----
write.xlsx2(tabulateSummary(f_fullVars),
                  "files/model outputs/tabulatedPGLS.xlsx", sheetName = "full")
write.xlsx2(tabulateSummary(f_redVars),
                  "files/model outputs/tabulatedPGLS.xlsx", sheetName = "red", append=T)
