## Spatial analysis, relating contribution of weather scores to geographic
## position

# housekeeping ----
library(dplyr); library(lme4); library(pbkrtest); library(ggplot2)

# read datasets ----
df_red <- readRDS("files/model outputs/redModels_explanatoryCapacity_spatial.rds")
df_full <- readRDS("files/model outputs/fullModels_explanatoryCapacity_spatial.rds")
df_both <- full_join(df_red, df_full, by=c("species", "gridID_200k"), suffix=c("_red", "_full")) %>%
    group_by(species) %>% mutate(meanImp=mean(modImp_red)) %>% ungroup %>%
    arrange(1-meanImp) %>%
    mutate(species=factor(species, levels=unique(species)))

suppPlot_1<- ggplot(df_both, aes(1-modImp_red, species)) + geom_point() +
    scale_x_continuous(breaks=seq(-1, 1, .05)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text.y  = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.x=element_text(hjust=70), 
          plot.title = element_text(face="bold", size=10)) +
    labs(title= "(b) full variable set", x=bquote("Contribution of weather,  "*1-MSE[w]^2/MSE[null]^2))

suppPlot_2<- ggplot(df_both, aes(1-modImp_full, species)) + geom_point() +
    scale_x_continuous(breaks=seq(-1, 1, .05)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size=7), 
          axis.title.y = element_blank(), 
          plot.title = element_text(face="bold", size=10)) +
    labs(title="(a) reduced variable set", x="") 

p_all <- egg::ggarrange(suppPlot_2, suppPlot_1, ncol=2)
ggsave("figures/testPlot_all.png", plot=p_all, height=297, width=190, units="mm") 
ggsave("figures/pdf/testPlot_all.pdf", plot=p_all, height=297, width=190, units="mm") 

## run model comparisons ----
# reduced variable set
lat_red_null <- lmer(1-modImp_red ~ poly(lon_red, 2) + (1|species) + (1|gridID_200k), df_both)
lat_red_1 <- lmer(1-modImp_red ~ lat_red + poly(lon_red, 2) + (1|species) + (1|gridID_200k), df_both)
lat_red_2 <- lmer(1-modImp_red ~ poly(lat_red, 2) + poly(lon_red, 2) + (1|species) + (1|gridID_200k), df_both)
KRmodcomp(lat_red_1, lat_red_null)
KRmodcomp(lat_red_2, lat_red_1)

lon_red_null <- lmer(1-modImp_red ~ poly(lat_red, 2) + (1|species) + (1|gridID_200k), df_both)
lon_red_1 <- lmer(1-modImp_red ~ lon_red + poly(lat_red, 2) + (1|species) + (1|gridID_200k), df_both)
lon_red_2 <- lmer(1-modImp_red ~ poly(lon_red, 2) + poly(lat_red, 2) + (1|species) + (1|gridID_200k), df_both)
KRmodcomp(lon_red_1, lon_red_null)
KRmodcomp(lon_red_2, lon_red_1)

# full variable set
lat_full_null <- lmer(1-modImp_full ~ poly(lon_red, 2) + (1|species) + (1|gridID_200k), df_both)
lat_full_1 <- lmer(1-modImp_full ~ lat_red + poly(lon_red, 2) + (1|species) + (1|gridID_200k), df_both)
lat_full_2 <- lmer(1-modImp_full ~ poly(lat_red, 2) + poly(lon_red, 2) + (1|species) + (1|gridID_200k), df_both)
KRmodcomp(lat_full_1, lat_full_null)
KRmodcomp(lat_full_2, lat_full_1)

lon_full_null <- lmer(1-modImp_full ~ poly(lat_red, 2) + (1|species) + (1|gridID_200k), df_both)
lon_full_1 <- lmer(1-modImp_full ~ lon_red + poly(lat_red, 2) + (1|species) + (1|gridID_200k), df_both)
lon_full_2 <- lmer(1-modImp_full ~ poly(lon_red, 2) + poly(lat_red, 2) + (1|species) + (1|gridID_200k), df_both)
KRmodcomp(lon_full_1, lon_full_null)
KRmodcomp(lon_full_2, lon_full_1)

# get model fits for significant models
preds_full <- df_both %>% na.omit %>%
    with(., data_frame(lon_red = seq(min(lon_red), max(lon_red, .5)), 
                       lat_red=median(lat_red))) %>%
    mutate(p = predict(lon_full_1, newdata=., re.form=NA))

pred_full <- data_frame(lon_red = seq(min(df_both$lat_red)))
preds_red <- df_both %>% 
    with(., data_frame(lon_red = seq(min(lon_red), max(lon_red, .5)), 
                       lat_red=median(lat_red))) %>%
             mutate(p = predict(lon_red_1, newdata=., re.form=NA))

# plot ----
p_red_lat <- ggplot(df_both, aes(lat_red, 1-modImp_red)) + geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=seq(-20, 70, 5)) + 
    scale_y_continuous(breaks=seq(-1, 1, .05)) + 
    labs(title = "(c) reduced variable set", 
         y= bquote("Contribution of weather,  "*1-MSE[w]^2/MSE[null]^2), x="Grid-cell latitude")  +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          plot.title=element_text(face="italic", size=11))

p_red_lon <- ggplot(df_both, aes(lon_red, 1-modImp_red)) + geom_point() +
    geom_line(data=preds_red, aes(y=p), size=1,  col="red") +
    theme_bw() +
    scale_x_continuous(breaks=seq(-20, 70, 5)) + 
    scale_y_continuous(breaks=seq(-1, 1, .05)) + 
    labs(title = "(d) reduced variable set", 
         y= bquote("Contribution of weather,  "*1-MSE[w]^2/MSE[null]^2), x="Grid-cell longitude")  +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          plot.title=element_text(face="italic", size=11))

preds_full <- df_both %>% na.omit %>%
    with(., data_frame(lon_red = seq(min(lon_red), max(lon_red, .5)), 
                       lat_red=median(lat_red))) %>%
    mutate(p = predict(lon_full_1, newdata=., re.form=NA))
max(preds_full$p) - min(preds_full$p)
max(preds_red$p) - min(preds_red$p)

p_full_lat <- ggplot(df_both %>% na.omit, aes(lat_red, 1-modImp_full)) + geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=seq(-20, 70, 5)) + 
    scale_y_continuous(breaks=seq(-1, 1, .05)) + 
    labs(title = "(a) full variable set", 
         y= bquote("Contribution of weather,  "*1-MSE[w]^2/MSE[null]^2), x="Grid-cell latitude")  +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          plot.title=element_text(face="italic", size=11))
p_full_lon <- ggplot(df_both %>% na.omit, aes(lon_red, 1-modImp_full)) + geom_point() +
    geom_line(data=preds_full, aes(y=p), size=1,  col="red") +
    scale_x_continuous(breaks=seq(-20, 70, 5)) + 
    scale_y_continuous(breaks=seq(-1, 1, .05)) + 
    theme_bw() +
    labs(title = "(b) full variable set", 
         y= bquote("Contribution of weather,  "*1-MSE[w]^2/MSE[null]^2), x="Grid-cell longitude")  +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour=NA),
          plot.title=element_text(face="italic", size=11))


p_all4 <- egg::ggarrange(p_full_lat, p_red_lat, p_full_lon, p_red_lon)         
ggsave("figures/suppPlots_latLon.png", plot=p_all4, width=200, height=200, units="mm")
ggsave("figures/pdf/suppPlots_latLon.pdf", plot=p_all4, width=200, height=200, units="mm")