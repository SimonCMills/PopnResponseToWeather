## read in and format trait DB for subsequent analysis. 
##
## Trait DB acquired from from: 
##  Storchová & Hořák (2018) Life-history characteristics of European birds. 
##  Global Ecology and Biogeography, 27, 400-406.

## dataset with all species in (note that the "full variable models" use a subset 
## of the species- as they have larger #dp require)
dataset <- readRDS("files/model outputs/redModels_explanatoryCapacity.rds")

speciesKey <- readRDS("files/birdDF.rds") %>%
    ungroup %>%
    dplyr::select(species, species_LB) %>%
    unique %>%
    left_join(dataset, ., by="species") %>%
    dplyr::select(species, species_LB)

dataset_wKey <- left_join(dataset, speciesKey)


traitDB <- read.table("files/Life-history characteristics of European birds.txt", "\t", header=T)
trait_cleaned <- traitDB %>%
    mutate(multibrood = ifelse(Broods.per.year==1, 0, 1), 
           LD_migrant = factor(Long.distance.migrant), 
           life_span = Life.span,
           hab_forest = ifelse(Deciduous.forest==1|Coniferous.forest==1|Woodland==1, 1, 0), 
           hab_open = ifelse(Shrub==1|Savanna==1|Tundra==1|Grassland==1|Mountain.meadows==1|Desert==1, 1, 0), 
           hab_aquatic = ifelse(Reed==1|Swamps==1|Freshwater==1|Marine==1, 1, 0), 
           diet_plant = ifelse(Folivore_Y==1|Frugivore_Y==1|Granivore_Y==1, 1,0), 
           diet_invert = ifelse(Arthropods_Y==1|Other.invertebrates_Y==1, 1, 0), 
           diet_vert = ifelse(Fish_Y==1|Other.vertebrates_Y==1, 1, 0), 
           diet_carrion = Carrion_Y) %>%
    dplyr::select(species_LB = Species, 
           body_mass = WeightU_MEAN, 
           multibrood, 
           LD_migrant,
           life_span,
           hab_forest, 
           hab_open, 
           hab_aquatic, 
           diet_plant, 
           diet_invert, 
           diet_vert, 
           diet_carrion)

## check correspondence between latin binomials (i.e. are species dropped b/c
## of scientific name mismatch) -- yes, 15 cases (and 4, for which #broods is 
## unknown or lifespan is unknown)
missingTrait <- left_join(speciesKey, trait_cleaned) %>% filter(!complete.cases(.)) %>%
    arrange(species_LB) %>%
    filter(!species %in% c("Cuckoo", "Fan-tailed Warbler", "Firecrest", "Little Bustard"))

(updatedLB <- missingTrait %>% 
    dplyr::select(species, species_LB) %>% 
    mutate(species_traitLB = c("Spatula clypeata", 
                          "Spatula querquedula", 
                          "Mareca strepera",
                          NA, 
                          "Linaria cannabina", 
                          "Acanthis flammea", 
                          "Spinus spinus", 
                          "Larus ridibundus", 
                          NA, 
                          "Leiopicus medius", 
                          "Dryobates minor", 
                          "Cyanecula svecica", 
                          "Poecile montanus", 
                          "Saxicola torquatus", 
                          "Lyrurus tetrix")))

trait_subset <- trait_cleaned %>%
    left_join(speciesKey, .) %>%
    filter(!species_LB %in% missingTrait$species_LB) %>%
    bind_rows(., updatedLB) %>%
    mutate(species_LB = ifelse(is.na(species_traitLB), species_LB, species_traitLB)) %>%
    dplyr::select(species, species_LB) %>%
    left_join(., trait_cleaned)

# now just 6 species that are not in phylogeny or miss trait info
trait_subset %>%
    filter(!complete.cases(.))

traitDB_final <- trait_subset %>%
    inner_join(., speciesKey, by="species", suffix=c("_trait", "")) %>%
    dplyr::select(species, species_LB, body_mass:diet_carrion) %>%
    filter(complete.cases(.))

## save ----
saveRDS(traitDB_final, "files/traitDB.rds")