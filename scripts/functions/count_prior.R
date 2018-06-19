# function to get count in year prior (checks that year interval ==1)
count_prior <- function(count, diff_priorYr, interval=1) {
    out <- rep(NA, length(count))
    # out vector with an observation in year prior (when interval == 1) gets the 
    # count from the year prior. When there is no observation in the year prior, 
    # it returns NA
    out[which(diff_priorYr==interval)] <- count[(which(diff_priorYr==interval)-interval)]
    out
}

# Check that this method performs correctly (TRUE)
# library(dplyr)
# check1 <- expand.grid(species = LETTERS[1:2], siteID = paste0("s", 1:2), year=2000:2006) %>%
#     as_data_frame %>%
#     arrange(species, siteID, year) %>%
#     mutate(count_t=year)
# 
# check1 %>%
#     filter(!year == 2002) %>%
#     group_by(species, siteID) %>%
#     mutate(diff_priorYr = year - c(NA, year[-length(year)])) %>%
#     View(.)
# 
# check1 %>%
#     filter(!year == 2002) %>%
#     group_by(species, siteID) %>%
#     mutate(diff_priorYr = year - c(NA, year[-length(year)]),
#            count_pre = count_prior(count_t, diff_priorYr, 1)) %>%
#     View(.)
