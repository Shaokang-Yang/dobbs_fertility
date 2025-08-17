library(tidyverse)
library(tidybayes)
library(posterior)
library(jsonlite)
library(kableExtra)
library(gt)

suffix <- "through_june"
fig_format <- "png"
## fig_prefix <- paste(suffix, "/", sep="")
fig_prefix <- ""

df <- read_csv(sprintf("/Users/shaokangyang/Library/CloudStorage/GoogleDrive-sky.ang510@gmail.com/My Drive/Code/fertility_results/2024/df_%s.csv", suffix))
df |>
  mutate(
    start_date = ym(paste(year, "-", bmcode * 2 - 1)),
    end_date = start_date + months(2) - days(1)
  ) -> df


fill_in_missing_denoms <- function(dat) {
  pop_index_2022 <- which.max(dat$year == 2022)
  pop_index_2021 <- which.max(dat$year == 2021)
  dat %>% mutate_at(vars(contains("pop")), ~ ifelse(is.na(.), .[pop_index_2022]^2 / .[pop_index_2021], .))
}

## Hacky imputation
df <- df %>%
  group_by(state) %>%
  group_modify(~ fill_in_missing_denoms(.)) %>%
  ungroup()

df$time <- lubridate::ymd(paste0(df$year, "-", 2*df$bmcode-1, "-01"))
df <- df %>% filter(time <= "2024-07-01" & time >= "2016-01-01") 
df$dobbs_code <- df$dobbscodev2
df <- df %>% group_by(state) %>% fill(exposed_births, .direction="down") %>% ungroup()
#df <- df %>% group_by(state) %>% fill(exposed_infdeaths, .direction="down") %>% ungroup()


relative_birth_rate <- df %>% group_by(state) %>% 
  mutate(ban = ifelse(any(exposed_births == 1), "Exposed (excl. Texas)", "Unexposed")) %>%
  mutate(ban = ifelse(is.na(ban), "Unexposed", ban)) %>% 
  mutate(ban = ifelse(state == "Texas", "Texas", ban)) %>%
  group_by(ban, time) %>% 
  summarize(births_total = sum(births_total), pop_total=sum(pop_total), exposed_births=mean(exposed_births)) %>% 
  ungroup() %>% group_by(ban) %>% 
  mutate(mean_br = mean(births_total[time < "2022-03-01"]/pop_total[time < "2022-03-01"])) %>%
  mutate(birthrate = (births_total/pop_total)/mean_br) 

ggplot(data=relative_birth_rate) + 
  geom_smooth(aes(x=time, y=birthrate, group=ban, col=ban), se=FALSE) + 
  geom_point(aes(x=time, y=birthrate, col=ban), alpha=0.5, data=relative_birth_rate) +
  theme_bw(base_size=16) + 
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        #legend.background = element_blank(),  # Make legend background transparent
        legend.title = element_blank()  ) +
  scale_color_manual(values=c("red", "orange", "dark gray"),
                     labels = c("Banned States (excl. Texas)", "Texas", "Non-banned States")) + 
  ylab("Relative Fertility Rate") + xlab("Year") + 
  geom_vline(xintercept=lubridate::date("2022-02-01"), color="orange", linetype="dashed") +
  geom_vline(xintercept=lubridate::date("2022-12-01"), color="red", linetype="dashed") 

ggsave(sprintf("figs/%smain_figures/relative_birthrate.%s", fig_prefix, fig_format), width=8, height=5)


df %>% group_by(state) %>%
  mutate(ban = ifelse(any(exposed_births == 1), "Exposed (excl. Texas)", "Unexposed")) %>%
  mutate(ban = ifelse(is.na(ban), "Unexposed", ban)) %>% 
  mutate(ban = ifelse(state == "Texas", "Texas", ban)) %>%
  group_by(ban, time) %>% 
  summarize(births_total = sum(births_total), pop_total=sum(pop_total), exposed_births=mean(exposed_births)) %>% 
  mutate(birthrate = births_total/pop_total*1000*6) %>%
  ungroup() %>% 
  ggplot() + 
  geom_smooth(aes(x=time, y=birthrate, col=ban), span=0.5, se=FALSE) + 
  geom_line(aes(x=time, y=birthrate, col=ban), alpha=0.5) +
  theme_bw(base_size=16) + 
  theme(legend.position = c(0.75, 0.99), legend.justification = c(1, 1),
        #legend.background = element_blank(),  # Make legend background transparent
        legend.title = element_blank()  ) +
  scale_color_manual(values=c("red", "orange", "dark gray"),
                     labels=c("States with bans (excl. Texas)", "Texas", "States without bans")) + 
  
  ylab("Fertility rate") + xlab("Year") +
  scale_x_date(
    date_breaks = "1 year",  # Set axis labels at yearly intervals
    date_labels = "%Y",
    limits=c(as.Date("2015-12-01"), as.Date("2023-08-01"))
  ) + guides(linetype="none") +
  geom_vline(xintercept=lubridate::date("2022-01-01"), color="orange", linetype="dashed") +
  geom_vline(xintercept=lubridate::date("2023-01-01"), color="red", linetype="dashed") +
  geom_text(x=lubridate::date("2022-02-01"), y=68, label="Texas exposure in effect", vjust=1, 
            angle=90, color="orange", size=3.5) + 
  geom_text(x=lubridate::date("2023-02-01"), y=68, label="Other exposures in effect", vjust=1, 
            angle=90, color="red", size = 3.5) + 
  # 12.5))
  coord_cartesian(xlim=c(as.Date("2015-12-01"), as.Date("2023-04-01")))

ggsave(sprintf("figs/%smain_figures/absolute_fertility.%s", fig_prefix, fig_format), width=8, height=5)  

df <- df %>% group_by(state) %>% 
  mutate(ban = ifelse(any(exposed_births == 1), TRUE, FALSE))

categories_list <- list(age = c("age1524", "age2534", "age3544"), 
                        edu = c("nohs", "hs", "somecoll", "coll"),
                        #edu = c("hs_less", "somecoll_more"),
                        insurance = c("medicaid", "nonmedicaid"),
                        marital = c("married", "unmarried"),
                        race = c("nhwhite", "nhblack", "hisp", "otherraceeth"),
                        total = c("total"))

###########################################################################
categories_list <- list(age = c("age1524", "age2534", "age3544"), 
                        edu = c("nohs", "hs", "somecoll", "coll"),
                        #edu = c("hs_less", "somecoll_more"),
                        insurance = c("medicaid", "nonmedicaid"),
                        marital = c("married", "unmarried"),
                        race = c("nhwhite", "nhblack", "hisp", "otherraceeth"),
                        total = c("total"))

## Load all datasets
file_dir <- "/Users/shaokangyang/Library/CloudStorage/GoogleDrive-sky.ang510@gmail.com/My Drive/Code/fertility_results"
types <- c("race", "total", "age", "edu", "insurance", "marital")

## These are the optimal ranks after PPCs
ranks <- c(10, 7, 10, 10, 11, 8)

all_samples <- tibble()
for(i in 1:length(types)) {
  type <- types[i]
  model_rank <- ranks[i]
  print(sprintf("%s %i", type, model_rank))
  df2 <- df
  if(type == "marital")
    df2 <- df2 %>% filter(state != "California")
  
  categories <- categories_list[[type]]
  type_samples <- read_csv(sprintf("%s/NB_births_%s_%i_%s.csv", file_dir, type, model_rank, suffix))
  
  categories <- categories_list[[type]]
  
  agg_category_name = "total"
  merged_df <- merge_draws_and_data(df2, type_samples,  categories=categories, agg_category_name=agg_category_name)
  df_ban_no_tx_treat_time <- merged_df %>% 
    filter(!(state %in% c("Texas", "Ban States")),  exposure_code==1) %>% 
    pull(time) %>% min()
  df_ban_no_tx <- merged_df %>% filter(state == "Ban States" ) %>% 
    mutate(exposure_code = ifelse(time >= df_ban_no_tx_treat_time, 1, 0))
  df_tx <- merged_df %>% filter(state == "Texas")
  df_ban_no_tx <- df_ban_no_tx %>% 
    mutate(state = "Ban States (excl. Texas)",
           ban = TRUE, 
           exposure_code = ifelse(time >= df_ban_no_tx_treat_time, 1, 0),
           ypred = df_ban_no_tx$ypred - df_tx$ypred,
           mu = log(exp(df_ban_no_tx$mu) - exp(df_tx$mu)),
           mu_treated = log(exp(df_ban_no_tx$mu_treated) - exp(df_tx$mu_treated)),
           pop = df_ban_no_tx$pop - df_tx$pop,
           births = df_ban_no_tx$births - df_tx$births,
           D = max(merged_df$D) + 1
    )
  merged_df <- merged_df %>% bind_rows(df_ban_no_tx)
  
  
  merged_df$type <- type
  merged_df$rank <- model_rank
  all_samples <- bind_rows(all_samples, merged_df)
}

all_samples <- all_samples %>% mutate(category = fct_recode(category,
                                                            "Hispanic" = "hisp",
                                                            "Non-Hispanic Black" = "nhblack",
                                                            "Non-Hispanic White" = "nhwhite",
                                                            "Other" = "otherraceeth",
                                                            "Total" = "total",
                                                            "15-24" = "age1524",
                                                            "25-34" = "age2534",
                                                            "35-44" = "age3544",
                                                            "Less than high school" = "nohs",
                                                            "High school diploma" = "hs",
                                                            "Some college" = "somecoll",
                                                            "College degree" = "coll",
                                                            "Medicaid" = "medicaid",
                                                            "Non-Medicaid" = "nonmedicaid",
                                                            "Married" = "married",
                                                            "Unmarried" = "unmarried"
)) %>% 
  mutate(category = fct_relevel(category,
                                "15-24", "25-34", "35-44",
                                "Less than high school", "High school diploma", "Some college", "College degree", 
                                "Medicaid", "Non-Medicaid",
                                "Married", "Unmarried",
                                "Hispanic", "Non-Hispanic Black", "Non-Hispanic White", "Other",
                                "Total"
  ))

quantiles_df <- all_samples %>% group_by(category, type, state, time) %>%
  summarize(ypred_mean=mean(ypred), 
            ypred_lower=quantile(ypred, 0.025), ypred_upper=quantile(ypred, 0.975), 
            births=mean(births), 
            exposure_code = first(exposure_code),
            ban = first(ban)) %>% ungroup()