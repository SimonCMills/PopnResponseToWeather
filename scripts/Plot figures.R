## Plot figures 

## housekeeping ----
library(dplyr);library(egg); library(stringr); library(viridis)

# create df to link numeric months to calendar months
month_vec <- month.abb[c(5:12, 1:4)]
varNames <- data_frame(var = c(paste0("precip_", 5:16), paste0("temp_", 5:16)), 
                       var_name = c(paste0("precip_", month_vec), paste0("temp_", month_vec))) %>%
    mutate(var_name = factor(var_name, levels=var_name))


## read in model outputs ----
# model explanatory capacity
modelPerf_full <- readRDS("files/model outputs/fullModels_explanatoryCapacity.rds") 
modelPerf_red <- readRDS("files/model outputs/redModels_explanatoryCapacity.rds") 

# bind into single df
modelPerf <- bind_rows(modelPerf_full, modelPerf_red, .id="id") %>%
    mutate(id = c("(a) Full variable set", "(b) Reduced variable set")[as.numeric(id)])

# get correlation coefficient between two measures
getCor <- left_join(modelPerf_full %>% dplyr::select(species, var_quotient),
          modelPerf_red %>% dplyr::select(species, var_quotient), by="species")
cor(getCor[,2:3], method = "spearm")

# model coefficients
coefs_full <- readRDS("files/model outputs/fullModels_coefs.rds") %>%
    filter(type %in% c("temp", "precip")) %>%
    mutate(var = str_extract(term, "p_[0-9]+"),
           var = gsub("p", "", var), 
           var = paste0(type, var)) %>%
    left_join(., varNames, by="var")

coefs_red <- readRDS("files/model outputs/redModels_coefs.rds") %>%
    filter(type %in% c("temp", "precip")) %>%
    mutate(var = str_extract(term, "p_....."),
           var = gsub("p", "", var), 
           var_name = paste0(type, var))

coefs <- bind_rows(coefs_full, coefs_red, .id="id") %>%
    mutate(id = c("(a) Full variable set", "(b) Reduced variable set")[as.numeric(id)])


## Plots ----
## Model improvement ----
# get mean
imp <- modelPerf %>% 
    group_by(id) %>%
    summarise(meanImp = mean(1-var_quotient), residSD = mean(sigma_full))

h1 <- modelPerf %>%
    filter(n_obs >= 300) %>%
    ggplot(aes(1-var_quotient, y = ..count../tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
    geom_histogram(boundary=0, binwidth=.005, fill="grey55") +
    theme_bw() +
    geom_vline(xintercept=0, lty=2) +
    labs(x= bquote("Contribution of weather,  "*1-sigma[w]^2/sigma[null]^2), y="Relative Frequency")  +
    # scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) +
    # guides(fill=F) +
    geom_vline(data=imp, aes(xintercept=meanImp), col="red") +
    facet_wrap(~id, nrow=2) +
    scale_x_continuous(breaks=seq(0, 1.05, .05)) +
    scale_y_continuous(breaks=seq(0, 1, .1))
ggsave("figures/ModelImprovementPlots.png", plot=h1, width=130, height=150, units="mm")

## Residual variation ----
h2 <- modelPerf %>%
    ggplot(aes(sigma_full, y = ..count../tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
    geom_histogram(boundary=0, binwidth=.025, fill="grey55") +
    theme_bw() +
    labs(x= bquote("Predictive uncertainty,  "*sigma[w]), y="Relative Frequency")  +
    # scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold"),#, 
          axis.text.x = element_text(lineheight=.001)) + #margin = margin(t = 0, r = 0, b = 0, l = 0), hjust=.5)) +
    guides(fill=F) +
    geom_vline(data=imp, aes(xintercept=residSD), col="red") +
    facet_wrap(~id, nrow=2) +
    scale_x_continuous(lim=c(0.15, 1.3), breaks=seq(.25, 1.25, .25), 
                       labels=c(expression(atop("0.25", "(0.61"*mu*", 1.63"*mu*")")), 
                                expression(atop("0.50", "(0.38"*mu*", 2.66"*mu*")")),
                                expression(atop("0.75", "(0.23"*mu*", 4.35"*mu*")")),
                                expression(atop("1.00", "(0.14"*mu*", 7.10"*mu*")")),
                                expression(atop("1.25", "(0.09"*mu*", 11.59"*mu*")"))
                                )) +
    scale_y_continuous(breaks=seq(0, 1, .05))
ggsave("figures/residVar.png", plot=h2, width=140, height=150, units="mm")

## Pval hist ----
h3 <- ggplot(coefs, aes(pVal, y = ..count../tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
    geom_histogram(boundary=0, binwidth=.01, colour="black", fill="grey80") +
    facet_wrap(~id, ncol=1) +
    scale_y_continuous(breaks=seq(0, 1, .02)) +
    theme_bw() +
    labs(x= "P-value", y="Relative Frequency")  +
    scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) + geom_hline(yintercept=.01, lty=2) 
ggsave("figures/pVals.png", plot=h3, width=130, height=150, units="mm")


## supp info plots ----
## Pvals by period ----
suppH3 <- coefs %>%
    mutate(var_name = gsub("_", ": ", var_name)) %>%
    filter(., id == "(a) Full variable set") %>%
    ggplot(aes(pVal, y = ..count../tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
    geom_histogram(boundary=0, binwidth=.05, colour="black", fill="grey80") +
    facet_wrap(~var_name, ncol=6) +
    theme_bw() +
    labs(x= "P-value", y="Relative Frequency")  +
    scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) + geom_hline(yintercept=.05, lty=2) 

suppH4 <- filter(coefs, id != "(a) Full variable set") %>%
    left_join(., 
              data_frame(var = c("_brMin", "_brMax", "_owMin", "_owMax"),
                         full_name = c(": br. season min.", ": br. season max.", ": overwinter min.", ": overwinter max."))) %>%
    mutate(full_name = paste0(type, full_name)) %>%
    ggplot(aes(pVal, y = ..count../tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
    geom_histogram(boundary=0, binwidth=.05, colour="black", fill="grey80") +
    facet_wrap(~full_name, ncol=4) +
    theme_bw() +
    labs(x= "P-value", y="Relative Frequency")  +
    scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) + geom_hline(yintercept=.05, lty=2) 

ggsave("figures/pVals_S1.png", plot=suppH3, width=297, height=190, units="mm")
ggsave("figures/pVals_S2.png", plot=suppH4, width=297, height=190, units="mm")

## pVals & mPerf by species, full mod ----
fullMods_pVals_bySp <- coefs_full %>% group_by(species) %>%
    summarise(propP = sum(pVal< .01)/n()) %>%
    arrange(propP) %>%
    mutate(species=factor(species, levels=species))
fullMods_mPerf_bySp <- modelPerf_full %>% 
    mutate(species=factor(species, levels=levels(fullMods_pVals_bySp$species)))

suppPlot_fullMPerf <- ggplot(fullMods_mPerf_bySp, aes(1-var_quotient, species)) + geom_point() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text.y  = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    labs(x=bquote("Contribution of weather,  "*sigma[w]^2/sigma[null]^2))

suppPlot_fullPVal <- ggplot(fullMods_pVals_bySp, aes(propP, species)) + geom_point() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size=7), 
          axis.title.y = element_text(hjust=1, angle=0)) +
    labs(y="(a)", x="Proportion of coefficient pvals < .01") +
    geom_vline(xintercept=0.01, lty=2)

suppPlot_full_both <- ggarrange(suppPlot_fullPVal, suppPlot_fullMPerf, ncol=2)
ggsave("figures/pVal&perf_bySp_S3.png", plot=suppPlot_full_both, height=297, width=190, units="mm") 

## pVals & mPerf by species, red mod ----
redMods_pVals_bySp <- coefs_red %>% group_by(species) %>%
    summarise(propP = sum(pVal < .01)/n()) %>%
    arrange(propP) %>%
    mutate(species=factor(species, levels=species))
redMods_mPerf_bySp <- modelPerf_red %>% 
    mutate(species=factor(species, levels=levels(redMods_pVals_bySp$species)))


suppPlot_redMPerf <- ggplot(redMods_mPerf_bySp, aes(1-var_quotient, species)) + geom_point() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text.y  = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    labs(x=bquote("Contribution of weather,  "*1-sigma[w]^2/sigma[null]^2))

suppPlot_redPVal <- ggplot(redMods_pVals_bySp, aes(propP, species)) + geom_point() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size=7), 
          axis.title.y = element_text(hjust=1, angle=0)) +
    labs(y="(b)", x="Proportion of coefficient pvals < .01") +
    geom_vline(xintercept=0.01, lty=2)

suppPlot_red_both <- ggarrange(suppPlot_redPVal, suppPlot_redMPerf, ncol=2)
ggsave("figures/pVal&perf_bySp_redSet_S4.png", plot=suppPlot_full_both, height=297, width=190, units="mm") 

## UK robustness check ----
# read in UK files and merge
modelPerf_UK_ECAD <- readRDS("files/model outputs/UKModels_ECAD_explanatoryCapacity.rds")
modelPerf_UK_UKCP <- readRDS("files/model outputs/UKModels_UKCP_explanatoryCapacity.rds")
modelPerf_UK <- merge(modelPerf_UK_ECAD, modelPerf_UK_UKCP, by=c("species", "n_grid", "n_obs", "n_scheme"),
                      suffixes=c("_ECAD", "_UKCP")) %>%
    as_data_frame()

coefs_UK_ECAD <- readRDS("files/model outputs/UKModels_ECAD_coefs.rds") %>%
    filter(type %in% c("temp", "precip")) %>%
    mutate(var = str_extract(term, "p_[0-9]+"),
           var = gsub("p", "", var), 
           var = paste0(type, var)) %>%
    left_join(., varNames, by="var") %>%
    dplyr::select(species, type, var_name, estimate, pVal)

coefs_UK_UKCP <- readRDS("files/model outputs/UKModels_UKCP_coefs.rds") %>%
    filter(type %in% c("temp", "precip")) %>%
    mutate(var = str_extract(term, "p_[0-9]+"),
           var = gsub("p", "", var), 
           var = paste0(type, var)) %>%
    left_join(., varNames, by="var") %>%
    dplyr::select(species, type, var_name, estimate, pVal)

coefs_UK <- merge(coefs_UK_ECAD, coefs_UK_UKCP, by=c("species", "type", "var_name"),
                      suffixes=c("_ECAD", "_UKCP")) %>%
    as_data_frame()

# get averages for plotting
imp_ECAD <- modelPerf_UK_ECAD %>%
    summarise(var_quotient=mean(var_quotient), sigma_full = mean(sigma_full))
imp_UKCP <- modelPerf_UK_UKCP %>%
    summarise(var_quotient=mean(var_quotient), sigma_full = mean(sigma_full))   
imp2 <-  bind_rows(imp_ECAD, imp_UKCP)

rob1 <- ggplot(modelPerf_UK, aes(1-var_quotient_ECAD, y = ..count../sum(..count..))) +
    geom_histogram(aes(1-var_quotient_UKCP), boundary=0, binwidth=.01, fill=cols_test[1]) +
    geom_histogram(boundary=0, binwidth=.01, alpha=.6, fill=cols_test[7]) +
    theme_bw() +
    geom_vline(xintercept=0, lty=2) +
    labs(x= bquote("Contribution of weather,  "*sigma[w]^2/sigma[null]^2), y="Relative Frequency")  +
    scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) +
    guides(fill=F) +
    geom_vline(data=imp2, aes(xintercept=1-var_quotient), col=c(cols_test[7], cols_test[1])) +
    # facet_wrap(~id, nrow=2, scale="free_y") +
    scale_x_continuous(breaks=seq(0, 1.05, .05))

rob2 <- modelPerf_UK %>%
    ggplot(aes(sigma_full_ECAD, y = ..count../sum(..count..))) + 
    geom_histogram(aes(sigma_full_UKCP), boundary=0, binwidth=.05, fill=cols_test[1]) +
    geom_histogram(boundary=0, binwidth=.05, alpha=.6, fill=cols_test[7]) +
    theme_bw() +
    labs(x= bquote("Predictive uncertainty,  "*sigma[w]), y="Relative Frequency")  +
    scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) +
    guides(fill=F) +
    geom_vline(data=imp2, aes(xintercept=sigma_full), col=c("indianred", "black")) 

rob3 <- ggplot(coefs_UK, aes(pVal_ECAD, y = ..count../sum(..count..))) +
    geom_histogram(aes(pVal_UKCP), boundary=0, binwidth=.01, fill=cols_test[1]) +
    geom_histogram(boundary=0, binwidth=.01, alpha=.6, fill=cols_test[7]) +
    # geom_histogram(aes(pVal_UKCP), boundary=0, binwidth=.01, alpha=.5) +
    theme_bw() +
    labs(x= "P-value", y="Relative Frequency")  +
    scale_fill_manual(values=c("grey80", "black")) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          strip.text=element_text(hjust=0, face="bold")) +
    geom_hline(yintercept=.01, lty=2) 

rob_all <- ggarrange(rob3, rob1, rob2, ncol=1)
ggsave("figures/UKfigs.png", plot=rob_all, width=130, height=225, units="mm")