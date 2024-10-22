library(tidybayes)
library(kableExtra)
library(ggrepel)
library(patchwork)
library(gt)


merge_draws_and_data <- function(dat, samples,
                                 categories = c("black", "other", "hisp", "white"), 
                                 outcome_string="births", denom_string="pop", 
                                 agg_category_name="total", verbose=FALSE) {
  
  
  ## When running model fit, categories are ordered alphabetically
  categories <- sort(categories)
  #categories <- c(categories, agg_category_name)
  
  if(outcome_string == "births") {
    dat <- dat %>% mutate(exposure_code = exposed_births)
  } else {
    dat <- dat %>% mutate(exposure_code = exposed_infdeaths)
  }

  dat <- dat %>%
    pivot_longer(
      cols = starts_with(paste0(outcome_string, "_")) | starts_with(paste0(denom_string, "_")),
      names_to = c(".value", "category"),
      names_pattern = "^([^_]*)_(.*)$"
    ) %>%
    select(state, exposure_code, ban, time, start_date, end_date, category, !!outcome_string, !!denom_string) %>%
    filter(category %in% categories)
  
  dat <- dat %>%
    mutate(state = factor(state)) %>%
    mutate(category = factor(category, levels = categories)) %>%
    mutate(state_code = as.numeric(state)) %>%
    mutate(category_code = as.numeric(factor(category, levels = categories))) %>% 
    mutate(time_code = as.numeric(factor(time))) %>%
    arrange(state, category, time)
  
  ban_group <- dat %>%
    filter(ban == 1) %>%
    group_by(time, start_date, end_date, category, category_code, time_code) %>%
    reframe(across(contains(outcome_string) | contains(denom_string), function(x) sum(x, na.rm = TRUE)),
      exposure_code = ifelse(any(exposure_code == 1), 1, 0)
    ) %>%
    ungroup() %>%
    mutate(state = "Ban States", ban = 1) %>%
    mutate(category_code = as.numeric(factor(category, levels = categories))) %>%
    mutate(time_code = as.numeric(factor(time))) %>%
    mutate(state_code = max(dat$state_code) + 1) %>%
    arrange(state, category, time)
    
  dat <- bind_rows(dat, ban_group)
  
  if (!any(categories == agg_category_name)) {
    dat_totals <- dat %>%
      group_by(state, exposure_code, ban, time, start_date, end_date, state_code, time_code) %>%
      summarize(
        across({{ outcome_string }}, sum, .names = "{outcome_string}"),
        across({{ denom_string }}, sum, .names = "{denom_string}")
      ) %>%
      mutate(category = agg_category_name, category_code = max(dat$category_code) + 1)

    dat <- bind_rows(dat, dat_totals)
  }
  ban_indices <- dat %>%
    filter(state != "Ban States") %>%
    select(ban, state) %>%
    distinct() %>%
    ungroup() %>%
    mutate(ban_index = row_number()) %>%
    filter(ban == TRUE) %>%
    pull(ban_index)

  draws_mat <- samples %>% as_draws_matrix()
  
  if(verbose) { print("Getting ypred draws...") }
  start <- Sys.time()
  
  if (any(grep("^te", colnames(draws_mat)))) {
    merged_draws <- draws_mat %>% spread_draws(ypred[K, D, N], mu[K, D, N], te[K, D, N])
    merged_draws <- merged_draws %>% mutate(mu_treated = mu + te)
  } else {
    merged_draws <- draws_mat %>% spread_draws(ypred[K, D, N], mu[K, D, N])
    merged_draws <- merged_draws %>% mutate(mu_treated = mu)
  }
  if(verbose) { print(Sys.time()-start) }

  if (!any(categories == agg_category_name)) {
    if (verbose) {
      print("Adding aggregated category to ypred...")
    }
    start <- Sys.time()
    merged_draws %>%
      ungroup() %>%
      group_by(D, N, .draw) %>%
      summarize(
        ypred = sum(ypred),
        mu = log(sum(exp(mu))),
        mu_treated = log(sum(exp(mu_treated)))
      ) %>% 
      mutate(K = max(dat$category_code)) -> merged_agg

    merged_draws <- bind_rows(merged_draws, merged_agg)
    if (verbose) {
      print(Sys.time() - start)
    }
  }
  
  if(verbose) { print("Adding Ban category to ypred...") }
  start <- Sys.time()
  merged_draws %>%
    ungroup() %>%
    filter(D %in% ban_indices) %>%
    group_by(K, N, .draw) %>%
    summarize(ypred = sum(ypred), 
              mu=log(sum(exp(mu))), 
              mu_treated = log(sum(exp(mu_treated)))) %>%
    mutate(D = length(unique(dat$state))) -> ban_preds
  merged_draws <- bind_rows(merged_draws, ban_preds)
  if(verbose) { print(Sys.time() - start) }
  
  if(verbose) { print("Merging data and draws...") }
  start <- Sys.time()
  merged_draws <- merged_draws %>%
    left_join(dat, by = c("N" = "time_code", "D" = "state_code", "K" = "category_code")) %>%
    ungroup()
  if(verbose) { print(Sys.time() - start) }

  merged_draws
}

make_all_te_plots <- function(merged_df, quantiles_df, state_name = "Texas", category="all", target="births") {
  p1 <- make_state_fit_plot(quantiles_df, state_name, category=category, target=target)
  p2 <- make_gap_plot(quantiles_df, state_name, target=target, category=category)
  p3 <- make_births_histogram(merged_df, state_name = state_name, category=category, target=target)

  (p1 + p2) / p3
}

make_state_fit_plot <- function(quantiles_df, state_name, category="total", target="births") {

  quantiles_df <- quantiles_df %>% filter(category == !!category)

  state_plot <- quantiles_df %>% filter(state == !!state_name) %>% ggplot() + 
  geom_point(aes(x=time, y=.data[[target]])) + 
  geom_ribbon(aes(x=time, ymin=ypred_lower, ymax=ypred_upper), alpha=0.5) + 
  geom_line(aes(x=time, y=ypred_mean), color="red") + theme_bw() + 
  xlab("Date") + ggtitle(state_name)
  if(target == "births") {
    state_plot <- state_plot + ylab("Births")
  } else {
    state_plot <- state_plot + ylab("Deaths")
  }
  
  

  if((quantiles_df$ban[quantiles_df$state == state_name])[1]){
    quantiles_df %>% filter(state == !!state_name, exposure_code == 1) %>% 
    summarize(treatment_date = first(time)) %>% 
    pull(treatment_date) -> treatment_date
    
    state_plot <- state_plot + geom_vline(xintercept=treatment_date, linetype="dashed")
  }
  
  state_plot
}

make_gap_plot <- function(quantiles_df, state_name = "Texas", category="all", target="births") {
  
  quantiles_df <- quantiles_df %>% filter(category == !!category & state == !!state_name)
  
  if ((quantiles_df$ban[quantiles_df$state == state_name])[1]) {
    quantiles_df %>%
      filter(state == !!state_name, exposure_code == 1) %>%
      summarize(treatment_date = first(time)) %>%
      pull(treatment_date) -> treatment_date
  }

  quantiles_df %>% 
  mutate(pre_ban = time < treatment_date) %>%
  ggplot() +
    geom_ribbon(aes(
      x = time,
      ymax = .data[[target]] / ypred_lower - 1,
      ymin = .data[[target]] / ypred_upper - 1,
      group = pre_ban
    ), alpha = 0.25) +
    geom_line(aes(x = time, y = .data[[target]] / ypred_mean - 1, group=pre_ban), color = "red") +
    theme_bw() +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed", alpha=0.75) +
    xlab("Date") -> gap_plot

  if(target == "births") {
    gap_plot <- gap_plot + ylab("Observed / Predicted Births - 1")
  } else {
    gap_plot <- gap_plot + ylab("Observed / Predicted Deaths - 1")
  }
  
  
  if((quantiles_df$ban[quantiles_df$state == state_name])[1]){
    gap_plot <- gap_plot + geom_vline(xintercept=treatment_date, linetype="dashed")
  }
  
  gap_plot
}

make_births_histogram <- function(merged_df, state_name = "Texas", category="all", 
                                  treatment_date = NULL, target="births") {
  
  merged_df <- merged_df %>%
    filter(
      category == !!category,
      state == !!state_name
  )

  if(is.null(treatment_date)) {
    merged_df %>% filter(exposure_code == 1) %>% 
    summarize(treatment_date = first(time)) %>% 
    pull(treatment_date) -> treatment_date
  }

  totals_df <- merged_df %>%
      filter(time >= treatment_date) %>%
      group_by(.draw) %>%
      summarize(obs = sum(.data[[target]]), pred = sum(ypred)) %>%
      mutate(diff = obs - pred)

  pval <- round(mean(totals_df$pred > totals_df$obs), 2)

  hist_plot <- totals_df %>% ggplot() + 
    geom_histogram(aes(x=diff), bins=50) + 
    geom_vline(xintercept=0, col="red", linetype="dashed") + 
    theme_bw() 
  if(target == "births") {
    hist_plot <- hist_plot + 
    ggtitle("Difference in Observed and Predicted Total Births") +
    xlab("Births")
  } else{
    hist_plot <- hist_plot + 
    ggtitle("Difference in Observed and Predicted Total Deaths") +
    xlab("Deaths")
  }
    
  hist_plot + annotate("text",
    x = Inf, y = Inf,
    label = sprintf("Pval = %.2f", pval),
    hjust = 1.2, vjust = 1.8
  )

}


#' Generate violin plots of births based on posterior predicted values
#'
#' This function generates violin plots of births based on predicted values.
#' It takes in a dataset, predicted draws, and optional parameters such as states, treatment date,
#' categories, target variable, and denominator variable.
#'
#' @param dat The dataset containing the birth data
#' @param pred_draws The predicted draws for births
#' @param states Optional vector of states to include in the plots
#' @param treatment_date Optional treatment date for filtering the data
#' @param categories Optional vector of categories to include in the plots
#' @param target The target variable (usually `births`)
#' @param denom The denominator variable for births
#'
#' @return A ggplot object containing the violin plots of births
#'
make_violins <- function(merged_df, states=NULL, treatment_date=NULL,
                         group_var = "state",
                         categories=NULL, target="births", denom="pop",
                         rate_normalizer=1000,
                         estimand = "diff", 
                         method="pred"){
  
  
  if(is.null(states)) { 
    states <- unique(merged_df$state[merged_df$ban == 1])
  }
  
  if (is.null(categories)) {
    categories <- unique(merged_df$category)
  } else {
    merged_df <- merged_df %>% filter(category %in% categories) 
  } 
  
  merged_df <- merged_df %>%
    filter(exposure_code == 1) %>%
    mutate(years = mean(interval(start_date, end_date) / years(1)))
  
  ## Remove Ban States, these will be recomputed below
  merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
  
  ## state %in% states, 
  
  if (method == "pred") {
    if (target == "births") {
      ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = mean(.data[[denom]])) %>%
        mutate(
          outcome_rate = (outcome / years) / (denom / rate_normalizer),
          ypred_rate = ypred / years / (denom / rate_normalizer)
        ) %>%
        ungroup() %>%
        group_by_at(c(".draw", group_var)) %>%
        summarize(
          causal_effect_diff = mean(outcome_rate - ypred_rate),
          causal_effect_ratio = mean(outcome_rate / ypred_rate),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>%
        ungroup()
    } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = sum(.data[[denom]])) %>%
        mutate(
          outcome_rate = outcome / (denom / rate_normalizer),
          ypred_rate = ypred / (denom / rate_normalizer)
        ) %>%
        ungroup() %>%
        group_by_at(c(".draw", group_var)) %>%
        mutate(
          causal_effect_ratio = mean(outcome_rate / ypred_rate),
          causal_effect_diff = mean(outcome_rate - ypred_rate),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>%
        ungroup()
    }
  } else {
    if (target == "births") {
      
      ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time", "state")) %>%
        mutate(treated = sum(exp(mu_treated)), untreated = sum(exp(mu)), denom = mean(.data[[denom]] * years)) %>%
        ungroup()
      
      ban_states_df <- state_df %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 1, state = "Ban States", ban = TRUE)
      ban_states_no_tx_df <- state_df %>%
        filter(state != "Texas") %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 2, state = "Ban States (excl. Texas)", ban = TRUE)
      
      
      state_df <- bind_rows(state_df %>% select_at(colnames(ban_states_df)), ban_states_df, ban_states_no_tx_df)
      state_df <- state_df %>% filter(state %in% states)
      # state_df <- bind_rows(state_df, ban_states_df, ban_states_no_tx_df)
      # state_df <- state_df %>% filter(state %in% states)
            
      state_df <- state_df %>% group_by_at(c(".draw", group_var)) %>%
        summarize(
          treated_rate = sum(treated) / sum(denom) * rate_normalizer,
          untreated_rate = sum(untreated) / sum(denom) * rate_normalizer,
          causal_effect_diff = treated_rate - untreated_rate,
          causal_effect_ratio = treated_rate / untreated_rate,
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>% ungroup()
        
    } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(treated = sum(exp(mu_treated)), untreated = sum(exp(mu)), denom = sum(.data[[denom]])) %>%
        ungroup()
        
      ban_states_df <- state_df %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 1, state = "Ban States", ban = TRUE)
      ban_states_no_tx_df <- state_df %>%
        filter(state != "Texas") %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 2, state = "Ban States (excl. Texas)", ban = TRUE)
      
      state_df <- bind_rows(state_df %>% select_at(colnames(ban_states_df)), ban_states_df, ban_states_no_tx_df)
      state_df <- state_df %>% filter(state %in% states)
      
        
      state_df <- state_df %>% group_by_at(c(".draw", group_var)) %>%
        summarize(
          treated_rate = sum(treated) / sum(denom) * rate_normalizer,
          untreated_rate = sum(untreated) / sum(denom) * rate_normalizer,
          causal_effect_diff = treated_rate - untreated_rate,
          causal_effect_ratio = treated_rate / untreated_rate,
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>% ungroup()
    }
  }

  state_df <- state_df %>%    
    mutate({{group_var}} := factor(.data[[group_var]])) %>%
    mutate({{group_var}} := fct_reorder(.data[[group_var]], causal_effect, .fun = median)) 

  stats_df <- state_df %>%
    group_by_at(group_var) %>%
    summarize(
      mean = round(mean(causal_effect), 3),
      pval = if(method=="pred"){ round(2*mean(untreated_rate > treated_rate), 3) } else { round(mean(treated_rate <= untreated_rate), 3) },
      Significance = factor(ifelse(pval < 0.05, "red", "black"), levels=c("red", "black"))
    ) %>%
    ungroup()
  


  vp <- state_df %>%
    ggplot(aes(x = .data[[group_var]], y = causal_effect)) +
    geom_violin(fill = "gray", alpha = 0.5, draw_quantiles = c(0.5)) +
    geom_hline(yintercept = ifelse(estimand=="diff", 0, 1), col = "red", linetype = "dashed") +
    geom_text(data = stats_df, aes(x = .data[[group_var]], y = Inf, label = pval, col = Significance), vjust = 2, fontface = "bold") +
    scale_colour_manual(values = c("red" = "red", "black" = "black")) +
    guides(colour = "none") +
    geom_text(data = stats_df, aes(x = .data[[group_var]], y = Inf, label = mean), vjust = 4, fontface = "bold") +
    theme_bw(base_size = 16) +    
    xlab(group_var) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (target == "births") {
    
    vp <- vp + labs(title = ifelse(estimand == "diff", "Post. Pred. Difference in birth rate", "Posterior Predictive Mult. change in birth rate"), subtitle = "Per 1,000 women per year")
    vp <- vp + ylab("Causal Effect")
  } else {
     vp <- vp + labs(title = ifelse(estimand == "diff", "Post. Pred. Difference in mortality rate", "Post. Pred. Mult. change in mortality rate"), subtitle = "Per 1,000 live births")
     vp <- vp + ylab("Mult. Change in deaths (per 1k live births)")
  }
  vp
}

#' Generate violin plots of births based on posterior predicted values
#'
#' This function generates violin plots of births based on predicted values.
#' It takes in a dataset, predicted draws, and optional parameters such as states, treatment date,
#' categories, target variable, and denominator variable.
#'
#' @param dat The dataset containing the birth data
#' @param pred_draws The predicted draws for births
#' @param states Optional vector of states to include in the plots
#' @param treatment_date Optional treatment date for filtering the data
#' @param categories Optional vector of categories to include in the plots
#' @param target The target variable (usually `births`)
#' @param denom The denominator variable for births
#'
#' @return A ggplot object containing the violin plots of births
#'
make_violins <- function(merged_df, states=NULL, treatment_date=NULL,
                         group_var = "state",
                         categories=NULL, target="births", denom="pop",
                         rate_normalizer=1000,
                         estimand = "diff", 
                         method="pred"){
  
  
  if(is.null(states)) { 
    states <- unique(merged_df$state[merged_df$ban == 1])
  }
  
  if (is.null(categories)) {
    categories <- unique(merged_df$category)
  } else {
    merged_df <- merged_df %>% filter(category %in% categories) 
  } 
  
  merged_df <- merged_df %>%
    filter(exposure_code == 1) %>%
    mutate(years = mean(interval(start_date, end_date) / years(1)))
  
  ## Remove Ban States, these will be recomputed below
  merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
  
  ## state %in% states, 
  
  if (method == "pred") {
    if (target == "births") {
      ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = mean(.data[[denom]])) %>%
        mutate(
          outcome_rate = (outcome / years) / (denom / rate_normalizer),
          ypred_rate = ypred / years / (denom / rate_normalizer)
        ) %>%
        ungroup() %>%
        group_by_at(c(".draw", group_var)) %>%
        summarize(
          causal_effect_diff = mean(outcome_rate - ypred_rate),
          causal_effect_ratio = mean(outcome_rate / ypred_rate),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>%
        ungroup()
    } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = sum(.data[[denom]])) %>%
        mutate(
          outcome_rate = outcome / (denom / rate_normalizer),
          ypred_rate = ypred / (denom / rate_normalizer)
        ) %>%
        ungroup() %>%
        group_by_at(c(".draw", group_var)) %>%
        mutate(
          causal_effect_ratio = mean(outcome_rate / ypred_rate),
          causal_effect_diff = mean(outcome_rate - ypred_rate),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>%
        ungroup()
    }
  } else {
    if (target == "births") {
      
      ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time", "state")) %>%
        mutate(treated = sum(exp(mu_treated)), untreated = sum(exp(mu)), denom = mean(.data[[denom]] * years)) %>%
        ungroup()
      
      ban_states_df <- state_df %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 1, state = "Ban States", ban = TRUE)
      ban_states_no_tx_df <- state_df %>%
        filter(state != "Texas") %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 2, state = "Ban States - No Texas", ban = TRUE)
      
      
      state_df <- bind_rows(state_df %>% select_at(colnames(ban_states_df)), ban_states_df, ban_states_no_tx_df)
      state_df <- state_df %>% filter(state %in% states)
      # state_df <- bind_rows(state_df, ban_states_df, ban_states_no_tx_df)
      # state_df <- state_df %>% filter(state %in% states)
            
      state_df <- state_df %>% group_by_at(c(".draw", group_var)) %>%
        summarize(
          treated_rate = sum(treated) / sum(denom) * rate_normalizer,
          untreated_rate = sum(untreated) / sum(denom) * rate_normalizer,
          causal_effect_diff = treated_rate - untreated_rate,
          causal_effect_ratio = treated_rate / untreated_rate,
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>% ungroup()
        
    } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(treated = sum(exp(mu_treated)), untreated = sum(exp(mu)), denom = sum(.data[[denom]])) %>%
        ungroup()
        
      ban_states_df <- state_df %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 1, state = "Ban States", ban = TRUE)
      ban_states_no_tx_df <- state_df %>%
        filter(state != "Texas") %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 2, state = "Ban States (excl. Texas)", ban = TRUE)
      
      state_df <- bind_rows(state_df %>% select_at(colnames(ban_states_df)), ban_states_df, ban_states_no_tx_df)
      state_df <- state_df %>% filter(state %in% states)
      
        
      state_df <- state_df %>% group_by_at(c(".draw", group_var)) %>%
        summarize(
          treated_rate = sum(treated) / sum(denom) * rate_normalizer,
          untreated_rate = sum(untreated) / sum(denom) * rate_normalizer,
          causal_effect_diff = treated_rate - untreated_rate,
          causal_effect_ratio = treated_rate / untreated_rate,
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>% ungroup()
    }
  }

  state_df <- state_df %>%    
    mutate({{group_var}} := factor(.data[[group_var]])) %>%
    mutate({{group_var}} := fct_reorder(.data[[group_var]], causal_effect, .fun = median)) 

  stats_df <- state_df %>%
    group_by_at(group_var) %>%
    summarize(
      mean = round(mean(causal_effect), 3),
      pval = if(method=="pred"){ round(2*mean(untreated_rate > treated_rate), 3) } else { round(mean(treated_rate <= untreated_rate), 3) },
      Significance = factor(ifelse(pval < 0.05, "red", "black"), levels=c("red", "black"))
    ) %>%
    ungroup()
  


  vp <- state_df %>%
    ggplot(aes(x = .data[[group_var]], y = causal_effect)) +
    geom_violin(fill = "gray", alpha = 0.5, draw_quantiles = c(0.5)) +
    geom_hline(yintercept = ifelse(estimand=="diff", 0, 1), col = "red", linetype = "dashed") +
    geom_text(data = stats_df, aes(x = .data[[group_var]], y = Inf, label = pval, col = Significance), vjust = 2, fontface = "bold") +
    scale_colour_manual(values = c("red" = "red", "black" = "black")) +
    guides(colour = "none") +
    geom_text(data = stats_df, aes(x = .data[[group_var]], y = Inf, label = mean), vjust = 4, fontface = "bold") +
    theme_bw(base_size = 16) +    
    xlab(group_var) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (target == "births") {
    
    vp <- vp + labs(title = ifelse(estimand == "diff", "Post. Pred. Difference in birth rate", "Posterior Predictive Mult. change in birth rate"), subtitle = "Per 1,000 women per year")
    vp <- vp + ylab("Causal Effect")
  } else {
     vp <- vp + labs(title = ifelse(estimand == "diff", "Post. Pred. Difference in mortality rate", "Post. Pred. Mult. change in mortality rate"), subtitle = "Per 1,000 live births")
     vp <- vp + ylab("Mult. Change in deaths (per 1k live births)")
  }
  vp
}

#' Generate violin plots of births based on posterior predicted values
#'
#' This function generates violin plots of births based on predicted values.
#' It takes in a dataset, predicted draws, and optional parameters such as states, treatment date,
#' categories, target variable, and denominator variable.
#'
#' @param dat The dataset containing the birth data
#' @param pred_draws The predicted draws for births
#' @param states Optional vector of states to include in the plots
#' @param treatment_date Optional treatment date for filtering the data
#' @param categories Optional vector of categories to include in the plots
#' @param target The target variable (usually `births`)
#' @param denom The denominator variable for births
#'
#' @return A ggplot object containing the violin plots of births
#'
make_interval_plot <- function(merged_df, states=NULL, treatment_date=NULL,
                         group_var = "state",
                         categories=NULL, target="births", denom="pop",
                         rate_normalizer=1000,
                         estimand = "diff", 
                         method="pred"){
  
  
  if(is.null(states)) { 
    states <- unique(merged_df$state[merged_df$ban == 1])
  }
  
  if (is.null(categories)) {
    categories <- unique(merged_df$category)
  } else {
    merged_df <- merged_df %>% filter(category %in% categories) 
  } 
  
  merged_df <- merged_df %>%
    filter(exposure_code == 1) %>%
    mutate(years = mean(interval(start_date, end_date) / years(1)))
  
  ## Remove Ban States, these will be recomputed below
  merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
  
  ## state %in% states, 
  
  if (method == "pred") {
    if (target == "births") {
      ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
      state_df <- merged_df %>%
        group_by_at(c(group_var, ".draw", "time")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = mean(.data[[denom]])) %>%
        mutate(
          outcome_rate = (outcome / years) / (denom / rate_normalizer),
          ypred_rate = ypred / years / (denom / rate_normalizer)
        ) %>%
        ungroup() %>%
        group_by_at(c(group_var, ".draw")) %>%
        summarize(
          causal_effect_diff = mean(outcome_rate - ypred_rate),
          causal_effect_ratio = mean(outcome_rate / ypred_rate),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>%
        ungroup()
    } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(group_var, ".draw", "time")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = sum(.data[[denom]])) %>%
        mutate(
          outcome_rate = outcome / (denom / rate_normalizer),
          ypred_rate = ypred / (denom / rate_normalizer)
        ) %>%
        ungroup() %>%
        group_by_at(c(".draw", group_var)) %>%
        mutate(
          causal_effect_ratio = mean(outcome_rate / ypred_rate),
          causal_effect_diff = mean(outcome_rate - ypred_rate),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>%
        ungroup()
    }
  } else {
    if (target == "births") {
      
      ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
      state_df <- merged_df %>%
        group_by_at(c(group_var, ".draw", "time", "state")) %>%
        mutate(treated = sum(exp(mu_treated)), untreated = sum(exp(mu)), denom = mean(.data[[denom]] * years)) %>%
        ungroup()
      
      ban_states_df <- state_df %>%
        group_by_at(c(".draw", setdiff(group_var, "state"), "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 1, state = "Ban States", ban = TRUE)
      ban_states_no_tx_df <- state_df %>%
        filter(state != "Texas") %>%
        group_by_at(c(".draw", setdiff(group_var, "state"), "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 2, state = "Ban States (excl. Texas)", ban = TRUE)
      
      
      state_df <- bind_rows(state_df %>% select_at(colnames(ban_states_df)), ban_states_df, ban_states_no_tx_df)
      state_df <- state_df %>% filter(state %in% states)
      # state_df <- bind_rows(state_df, ban_states_df, ban_states_no_tx_df)
      # state_df <- state_df %>% filter(state %in% states)
            
      state_df <- state_df %>% group_by_at(c(".draw", group_var)) %>%
        summarize(
          treated_rate = sum(treated) / sum(denom) * rate_normalizer,
          untreated_rate = sum(untreated) / sum(denom) * rate_normalizer,
          causal_effect_diff = treated_rate - untreated_rate,
          causal_effect_ratio = 100*(treated_rate / untreated_rate - 1),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>% ungroup()
        
    } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(".draw", group_var, "time")) %>%
        mutate(treated = sum(exp(mu_treated)), untreated = sum(exp(mu)), denom = sum(.data[[denom]])) %>%
        ungroup()
        
      ban_states_df <- state_df %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 1, state = "Ban States", ban = TRUE)
      ban_states_no_tx_df <- state_df %>%
        filter(state != "Texas") %>%
        group_by_at(c(".draw", "category", "time")) %>%
        summarize(treated = sum(treated), untreated = sum(untreated), denom = sum(denom)) %>%
        mutate(D = max(state_df$D) + 2, state = "Ban States (excl. Texas)", ban = TRUE)
      
      state_df <- bind_rows(state_df %>% select_at(colnames(ban_states_df)), ban_states_df, ban_states_no_tx_df)
      state_df <- state_df %>% filter(state %in% states)
      
        
      state_df <- state_df %>% group_by_at(c(".draw", group_var)) %>%
        summarize(
          treated_rate = sum(treated) / sum(denom) * rate_normalizer,
          untreated_rate = sum(untreated) / sum(denom) * rate_normalizer,
          causal_effect_diff = treated_rate - untreated_rate,
          causal_effect_ratio = 100*(treated_rate / untreated_rate - 1),
          causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)
        ) %>% ungroup()
    }
  }

  
  state_df <- state_df %>%
    mutate(state = factor(state)) %>%
    mutate(state = fct_relevel(state, "Ban States", "Ban States (excl. Texas)")) 

  stats_df <- state_df %>%
    group_by_at(group_var) %>%
    summarize(
      mean = round(mean(causal_effect), 3),
      pval = if(method=="pred"){ round(2*mean(untreated_rate > treated_rate), 3) } else { round(mean(treated_rate <= untreated_rate), 3) },
      Significance = factor(ifelse(pval < 0.05, "red", "black"), levels=c("red", "black"))
    ) %>%
    ungroup()
  
  
  color_group <- setdiff(group_var, "state")
  state_df %>%
    #filter(category != "total") %>%
    ggplot(aes(x = state, y = causal_effect, color = fct_rev(!!sym(color_group)))) +
    ggdist::stat_pointinterval(aes(
      alpha = after_stat(level)
    ), position = "dodge", .width = c(0.95, 0.67)) +
    ggdist::scale_interval_alpha_continuous(range = c(0.75, 1)) +
    colorspace::scale_color_discrete_qualitative() +
    scale_alpha_manual(values = c(0.5, 1), ) +
    theme_bw(base_size = 16) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Percent Change From Expected") +
    xlab("") +
    guides(colour = guide_legend(reverse=T)) + 
    coord_flip() + facet_grid(state ~ ., scales = "free_y") + theme(strip.text.y = element_blank())
    



}



make_violin_diffs <- function(merged_df, state="Texas", treatment_date=NULL, 
                                target="births", denom="pop",
                                rate_normalizer=1000,
                                estimand = "diff"){
    
  merged_df <- merged_df %>%
    filter(state == !!state, exposure_code == 1) %>%
    group_by(state) %>%
    mutate(years = interval(min(start_date, na.rm=TRUE), max(end_date, na.rm=TRUE)) / years(1)) %>%
    ungroup()
      
  if(target == "births") {
    ## Compute ratio of birth rates (rates per 1000 people per year (or the rate_normalizer))
    state_df <- merged_df %>%
      group_by_at(c(".draw", "category")) %>%
      mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = mean(.data[[denom]])) %>%
      mutate(
        outcome_rate = (outcome / years) / (denom / rate_normalizer),
        ypred_rate = ypred / years / (denom / rate_normalizer)
      ) %>%
      mutate(
        causal_effect_diff = outcome_rate - ypred_rate,
        causal_effect_ratio = outcome_rate / ypred_rate,
        causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)) %>%
      ungroup()
  } else {
      ## Compute difference in death rate per 1000 births (or the rate_normalizer)
      state_df <- merged_df %>%
        group_by_at(c(".draw", "category")) %>%
        mutate(outcome = sum(.data[[target]]), ypred = sum(ypred), denom = sum(.data[[denom]])) %>%
        mutate(outcome_rate = outcome / (denom / rate_normalizer),
               ypred_rate = ypred / (denom / rate_normalizer)) %>%
      mutate(causal_effect_ratio = outcome_rate / ypred_rate,
             causal_effect_diff = outcome_rate - ypred_rate,
             causal_effect = ifelse(estimand == "diff", causal_effect_diff, causal_effect_ratio)) %>%
      ungroup()
  }
  

  state_df <- state_df %>%
    mutate(category = factor(category)) %>%
    mutate(category = fct_reorder(category, causal_effect, .fun = median))

  cat_levels <- setdiff(levels(state_df$category), "total")
  diff_cats <- combn(cat_levels, 2) %>% t
  diff_cat1 <- diff_cats[, 2]
  diff_cat2 <- diff_cats[, 1]
   
  #  state_df <- state_df %>%    
  #   mutate({{group_var}} := factor(.data[[group_var]])) %>%
  #   mutate({{group_var}} := fct_reorder(.data[[group_var]], causal_effect, .fun = median)) 

  state_wide_df <- state_df %>%
    pivot_wider(
      id_cols = c(.draw, time), names_from = category,
      values_from = causal_effect
    )
  
  # state_df <- state_df %>%    
  #   mutate({{group_var}} := factor(.data[[group_var]])) %>%
  #   mutate({{group_var}} := fct_reorder(.data[[group_var]], causal_effect, .fun = median)) 

  for(i in 1:length(diff_cat1)) {
    state_wide_df <- state_wide_df %>% 
      mutate(!!paste0(diff_cat1[i], " - ", diff_cat2[i]) := !!sym(diff_cat1[i]) - !!sym(diff_cat2[i]))
  }
  
  
  state_long_df <- state_wide_df %>% pivot_longer(
    cols = contains(" - "),
    names_to = "diff_cat", values_to = "diff"
  )

  state_long_df <- state_long_df %>%
     mutate(diff_cat = factor(diff_cat)) %>%
     mutate(diff_cat := fct_reorder(diff_cat, diff, .fun = median))   

  stats_df <- state_long_df %>%
    group_by(diff_cat) %>%
    summarize(
      mean = round(mean(diff), 3),
      pval = round(mean(diff < 0), 3),
      Significance = factor(ifelse(pval < 0.05, "red", "black"), levels=c("red", "black"))
    ) %>%
    ungroup() 
  
  

  state_long_df %>%
    ggplot(aes(x = diff_cat, y = diff)) +
    geom_violin(fill = "gray", alpha = 0.5, draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    geom_text(data = stats_df, aes(x = diff_cat, y = Inf, label = pval, col = Significance), vjust = 2, fontface = "bold") +
      scale_colour_manual(values=c("red"="red", "black"="black")) + guides(colour="none") +
    geom_text(data=stats_df, aes(x=diff_cat, y=Inf, label=mean), vjust=4, fontface='bold') +
    theme_bw(base_size=16) + ggtitle("Difference in Group Effects") +
    xlab("State") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("") 


}

make_table <- function(merged_df, 
                       target_state = "Texas", target="births", denom="pop",
                       rate_normalizer=1000, plot_type="exploratory") {
  
  if(target_state == "Ban States") {
    merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
    merged_df <- merged_df %>%
      filter(exposure_code == 1) %>%
      ## Aggregate over all banned states
      group_by(category, .draw, time) %>% 
      summarise({{target}} := sum(.data[[target]]), 
                denom = sum(.data[[denom]]), 
                ypred=sum(ypred), 
                mu = log(sum(exp(mu))),
                mu_treated = log(sum(exp(mu_treated))),
                years=mean(interval(start_date, end_date) / years(1)))
      
  } else if(target_state == "Ban States (excl.Texas)") {
    merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
    merged_df <- merged_df %>%
      filter(state != "Texas") %>%
      filter(exposure_code == 1) %>%
      ## Aggregate over all banned states
      group_by(category, .draw, time) %>% 
      summarise({{target}} := sum(.data[[target]]), 
                denom = sum(.data[[denom]]), 
                ypred=sum(ypred), 
                mu = log(sum(exp(mu))),
                mu_treated = log(sum(exp(mu_treated))),
                years=mean(interval(start_date, end_date) / years(1)))

  } else {
    merged_df <- merged_df %>%
      filter(state == target_state, exposure_code == 1) %>%
      mutate(years = interval(start_date, end_date) / years(1), denom=.data[[denom]])
  }
  
  table_df <- merged_df %>%
    ungroup() %>%
    ## Aggregate over time
    group_by(category, .draw) %>%
    summarize(
      ypred = sum(ypred),
      outcome = sum(.data[[target]]), years = mean(years),
      treated = sum(exp(mu_treated)), untreated = sum(exp(mu)),
      denom = ifelse(target == "births", sum(denom * years, na.rm = TRUE), sum(denom, na.rm = TRUE)),
      treated_rate = treated / denom * rate_normalizer,
      untreated_rate = untreated / denom * rate_normalizer,
      outcome_rate = round(outcome / denom * rate_normalizer, 2),
      outcome_diff = round(treated - untreated)
    ) %>%
    ungroup() %>%
    ## Compute quantiles of effects
    group_by(category) %>%
    summarize(
      ypred_mean = mean(ypred),
      outcome = mean(outcome),
      outcome_diff_mean = round(mean(outcome_diff)), 
      outcome_diff_lower = round(quantile(outcome_diff, 0.025)), 
      outcome_diff_upper = round(quantile(outcome_diff, 0.975)),
      outcome_rate = mean(outcome_rate),
      ypred_lower = quantile(ypred, 0.025), ypred_upper = quantile(ypred, 0.975),
      treated_mean = mean(treated), treated_lower = quantile(treated, 0.025), treated_upper = quantile(treated, 0.975),
      untreated_mean = mean(untreated), untreated_lower = quantile(untreated, 0.025), untreated_upper = quantile(untreated, 0.975),      
      treated_rate_mean = mean(treated_rate), treated_rate_lower = quantile(treated_rate, 0.025), treated_rate_upper = quantile(treated_rate, 0.975),
      untreated_rate_mean = mean(untreated_rate), untreated_rate_lower = quantile(untreated_rate, 0.025), untreated_rate_upper = quantile(untreated_rate, 0.975), 
      causal_effect_diff_mean = mean(treated_rate - untreated_rate), causal_effect_diff_lower = quantile(treated_rate - untreated_rate, 0.025), causal_effect_diff_upper = quantile(treated_rate - untreated_rate, 0.975),
      causal_effect_ratio_mean = mean(treated_rate / untreated_rate), causal_effect_ratio_lower = quantile(treated_rate / untreated_rate, 0.025), causal_effect_ratio_upper = quantile(treated_rate / untreated_rate, 0.975),
      denom = mean(denom),
      pval = 2*mean(untreated > treated)
     )
    
    
  table_df <- table_df %>%
  mutate(
    # ypred_mean_rate = ypred_mean / years / (denom / rate_normalizer),
    rate_diff = round(causal_effect_diff_mean, 2),
    rate_diff_lower = round(causal_effect_diff_lower, 2),
    rate_diff_upper = round(causal_effect_diff_lower, 2),
    mult_change = round(causal_effect_ratio_mean, 3),
    mult_change_lower = round(causal_effect_ratio_lower, 3),
    mult_change_upper = round(causal_effect_ratio_upper, 3)
  )
  
  table_df <- table_df %>%        
    select(category, outcome, outcome_diff_mean, outcome_diff_lower, outcome_diff_upper, outcome_rate, 
           mult_change, rate_diff, mult_change_lower, mult_change_upper, pval) %>%
    relocate(c("pval"), .after = "category") %>%
    arrange(desc(mult_change))
  
  table_out <- table_df %>%
  gt() |>
  tab_spanner(
    label = "Count",
    columns = c("outcome", "outcome_diff_mean", "outcome_diff_lower", "outcome_diff_upper")
  ) |>
  tab_spanner(
    label = "Rate",
    columns = c("outcome_rate", "rate_diff", "mult_change", "mult_change_lower", "mult_change_upper")
  )
  if(target == "births") {
    table_out <- table_out %>% cols_label(
      category = "Category", outcome = "Births", 
      outcome_diff_lower = "Lower", outcome_diff_upper = "Upper", 
      outcome_diff_mean = "Birth Diff.", outcome_rate = "Birth Rate", 
      mult_change_lower = "Lower", mult_change_upper = "Upper", 
      rate_diff = "Rate Diff.",
      mult_change = "Rate Mult.",
      pval = "Pval"
    )
  } else {
    table_out <- table_out %>% cols_label(
      category = "Category", outcome = "Deaths", outcome_diff_lower = "Lower", 
      outcome_diff_upper = "Upper", outcome_diff = "Death Diff.",
      outcome_rate = "Death Rate", mult_change_lower = "Lower", 
      mult_change_upper = "Upper", rate_diff = "Rate Diff.",
      mult_change = "Rate Mult.",
      pval = "Pval"
    )
  }

  if (any(table_df$mult_change > 1)) {
   table_out <- table_out |> data_color(
      columns = mult_change,
      rows = mult_change > 1,
      method = "numeric",
      domain = c(1, max(table_df$mult_change)),
      palette = "OrRd"
    )
  }
  if(any(table_df$pval < 0.05)) {
    table_out <- table_out |>
      data_color(
        columns = pval,
        target_columns = c("category", "pval"),
        rows = pval < 0.05,
        method = "numeric",
        domain = c(0, .05),
        palette = "OrRd",
        reverse = TRUE
      )
  }
  table_out
}


make_all_states_plot <- function(dat, focus_state = "Texas") {

  focus_dat <- dat %>% filter(state == !!focus_state) %>% 
    mutate(Ef = Ef_samp, diff=births-Ef_samp, 
           Ef_lower=Ef_samp_lower, Ef_upper=Ef_samp_upper)
  
  gp <- dat %>%
    mutate(Ef = Ef_samp, diff=births-Ef_samp, 
           Ef_lower=Ef_samp_lower, Ef_upper=Ef_samp_upper) %>%
    group_by(state) %>% 
    mutate(mean_births = mean(births), mean_population = mean(population)) %>% ungroup() %>%
    ggplot(aes(x=year_mth, y=births)) + 
    geom_line(aes(y=diff/sqrt(population/1000), group=state), color="dark gray", alpha=0.5) +
    geom_hline(yintercept=0, color='black', linetype="dashed") +
    geom_vline(xintercept=2022.25, color='black', linetype="dashed") +
    labs(x="Date", y="Scaled Residual") + theme_bw()  +
    geom_line(data=focus_dat, 
              aes(y=diff/sqrt(population/1000), x=year_mth, 
                  color=ifelse(year_mth < 2022.25, "pre-treatment", "post-treatment"))) + 
    scale_color_manual(name="Texas", values=set1[1:2]) + 
    ggtitle("Residual Plot")

  gp

}



### Posterior Predictive Checks

make_acf_ppc_plot <- function(merged_df, lag=6, 
                              outcome="births",
                              categories = NULL) {
  
  if (is.null(categories)) {
    categories <- unique(merged_df$category)
  }

  ban_states <- merged_df %>%
    filter(ban == 1) %>%
    pull(state) %>%
    unique()

  acf_stats <- merged_df %>%
    filter(exposure_code == 0, state %in% ban_states) %>%
    filter(state != "Ban States") %>%
    filter(category %in% categories) %>%
    mutate(pred_diff = ypred - exp(mu)) %>%
    mutate(obs_diff = .data[[outcome]] - exp(mu)) %>%
    group_by(state, category, .draw) %>%
    summarise(
      obs_ac = acf(obs_diff, lag.max = lag, plot = FALSE)$acf[lag+1, 1, 1],
      pred_ac = acf(pred_diff, lag.max = lag, plot = FALSE)$acf[lag+1, 1, 1],
      diff_in_ac = obs_ac - pred_ac)

  pvals <- acf_stats %>%
      group_by(state, category) %>%
      summarize(pval = mean(diff_in_ac < 0)) %>%
      ungroup() %>%
      filter(category %in% categories) %>%
      filter(state %in% ban_states)

  acf_plt <- acf_stats %>%
    filter(category %in% categories) %>%
    ggplot() +
    geom_histogram(aes(x = diff_in_ac), alpha = 0.5) +
    geom_text(data=pvals, aes(label = round(pval, 3)), 
                y = Inf, x = Inf, hjust = 1, vjust = 1, col="red") +
    geom_vline(xintercept=0, col="red", linetype="dashed") + 
    ggtitle(sprintf("Difference in Residual Autocorrelation (Lag %i)", lag)) +
    facet_wrap(~state + category, scales="free", ncol=3) + theme_bw() + xlab("Observed - Predicted Autocorrelation")
  
  list("pvals"=pvals$pval, "acf_plt"=acf_plt)

} 

make_rmse_ppc_plot <- function(merged_df, 
                      outcome="births", categories = NULL) {
  

  if(is.null(categories)) {
    categories <- unique(merged_df$category)
  }
  rmse_stats <- merged_df %>%
    filter(exposure_code == 0) %>% 
    mutate(pred_diff = ypred - exp(mu)) %>%
    mutate(obs_diff = .data[[outcome]] - exp(mu)) %>%
    group_by(state, category, .draw) %>%
    summarise(
      rmse_pred_diff = sqrt(mean(pred_diff^2)),
      rmse_obs_diff = sqrt(mean(obs_diff^2))
    ) %>% mutate(diff_in_diff = rmse_obs_diff - rmse_pred_diff)
  
    ban_states <- merged_df %>%
      filter(ban == 1) %>%
      pull(state) %>%
      unique()

    pvals <- rmse_stats %>%
      group_by(state, category) %>%
      summarize(pval = mean(diff_in_diff < 0)) %>%
      ungroup() %>%
      filter(category %in% categories) %>%
      filter(state %in% ban_states)
    
    rmse_plt <- rmse_stats %>% 
      filter(state %in% ban_states) %>%
      filter(category %in% categories) %>%
      ggplot() +
      geom_histogram(aes(x = diff_in_diff), alpha = 0.5) +
      geom_text(data=pvals, aes(label = round(pval, 3)), 
                y = Inf, x = Inf, hjust = 1, vjust = 1, col="red") +
      geom_vline(xintercept = 0, col = "red", linetype = "dashed") +
      facet_wrap(~ state + category, scales="free", ncol=3) +
      theme_bw() +
      ggtitle("Difference in RMSE") + 
      xlab("Observed - Predicted RMSE")

    list("pvals"=pvals$pval, "rmse_plt"=rmse_plt)

}

## Maximum absolute residual

make_abs_res_ppc_plot <- function(merged_df,  
                                  outcome="births", categories = NULL) {
  
  if (is.null(categories)) {
    categories <- unique(merged_df$category)
  }
  max_stats <- merged_df %>%
    filter(exposure_code == 0) %>% 
    mutate(pred_diff = ypred - exp(mu)) %>%
    mutate(obs_diff = .data[[outcome]] - exp(mu)) %>%
    group_by(state, category, .draw) %>%
    summarise(
      max_pred_diff = max(abs(pred_diff)),
      max_obs_diff = max(abs(obs_diff))
    ) %>% mutate(diff_in_diff = max_obs_diff - max_pred_diff)
  
  ban_states <- merged_df %>%
    filter(ban == 1) %>%
    pull(state) %>%
    unique()

  pvals <- max_stats %>%
    group_by(state, category) %>%
    summarize(pval = mean(diff_in_diff < 0)) %>%
    ungroup() %>%
    filter(category %in% categories) %>%
    filter(state %in% ban_states)

  max_plt <- max_stats %>%
    filter(state %in% ban_states) %>%
    filter(category %in% categories) %>%
    ggplot() +
    geom_histogram(aes(x = diff_in_diff), alpha = 0.5) +
    geom_text(data=pvals, aes(label = round(pval, 3)), 
              y = Inf, x = Inf, hjust = 1, vjust = 1, col="red") +
    geom_vline(xintercept = 0, col = "red", linetype = "dashed") +
    facet_wrap(~ state + category, scales="free", ncol=3) +
    theme_bw() +
    ggtitle("Difference in Maximum Absolute Predicted Residual") + xlab("Observed - Predicted Max Residual")
  
  list("pvals"=pvals$pval, "max_plt"=max_plt)

}

make_unit_corr_ppc_plot <- function(merged_df,
                                    max_treat_date = "2022-04-01", 
                                    categories = NULL,
                                    ndraws_to_use=1000, outcome="births") {
  
  if(is.null(categories)) {
    categories <- unique(merged_df$category)
  }                           

  eval_stats <- merged_df %>%
    filter(time < max_treat_date) %>%
    filter(.draw < ndraws_to_use) %>%
    filter(category %in% categories) %>%
    mutate(obs_residual = .data[[outcome]] - exp(mu), pred_residual = ypred - exp(mu)) %>%
    group_by(state, category) %>%
    mutate(na_outcomes = mean(is.na(.data[[outcome]]))) %>%
    ungroup() %>%
    filter(na_outcomes < 0.25) %>%
    group_by(category, .draw) %>%
    summarise(
      obs_sval = sqrt(eigen(cor(matrix(obs_residual, ncol = length(unique(D))), use = "pairwise.complete.obs"))$values[1]),
        pred_sval = sqrt(eigen(cor(matrix(as.logical(obs_residual)*pred_residual, ncol = length(unique(D))), use = "pairwise.complete.obs"))$values[1])
      ) %>%
    mutate(eval_diff = obs_sval - pred_sval)
    
    pvals <- eval_stats %>%
      group_by(category) %>%
      summarize(pval = mean(eval_diff < 0)) %>%
      ungroup() %>%
      filter(category %in% categories)
    
    eval_plt <- eval_stats %>%
      ggplot() +
      geom_histogram(aes(x = eval_diff), alpha = 0.5) +
      geom_text(data=pvals, aes(label = round(pval, 3)), 
              y = Inf, x = Inf, hjust = 1, vjust = 1, col="red") +
      geom_vline(xintercept = 0, col = "red", linetype = "dashed") +
      facet_wrap(~category, scales = "free", ncol=2) +
      theme_bw() +
      ggtitle("Differnece in State Correlations") +
      xlab("Observed - Predicted Spectral Norm")
    

    list("pvals"=pvals$pval, "eval_plt"=eval_plt)
    # tibble(eval=evals_rep) %>% ggplot(aes(x=evals_rep)) + geom_histogram() + 
    #   geom_vline(xintercept=eval_obs, col="red", linetype="dashed") + 
    #    + theme_bw() +
    #   annotation_custom(grid::textGrob(sprintf("p(T > obs) = %s", format(round(mean(evals_rep > eval_obs), 3), nsmall=3)), 0.15, 0.95, 
    #                                    gp = grid::gpar(col = "red", fontsize = 20, hjust=1)))
}



make_fertility_table <- function(merged_df, 
                       target_state = "Texas", target="births", denom="pop",
                       rate_normalizer=1000, plot_type="exploratory",
                       tab_caption = "Table 1. Estimated difference in cumulative observed vs expected births (count and rate) in all states that banned abortion in months affected by bans (January 2023 through December 2023), overall and by socioeconomic characteristics") {

 if(target_state == "Ban States") {
    merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
    merged_df <- merged_df %>%
      filter(exposure_code == 1) %>%
      ## Aggregate over all banned states
      group_by(type, category, .draw, time) %>% 
      summarise({{target}} := sum(.data[[target]]), 
                denom = sum(.data[[denom]]), 
                ypred=sum(ypred), 
                mu = log(sum(exp(mu))),
                mu_treated = log(sum(exp(mu_treated))),
                years=mean(interval(start_date, end_date) / years(1)))
      
  } else if(target_state == "Ban States (excl.Texas)") {
    merged_df <- merged_df %>% filter(!state %in% c("Ban States", "Ban States (excl. Texas)"))
    merged_df <- merged_df %>%
      filter(state != "Texas") %>%
      filter(exposure_code == 1) %>%
      ## Aggregate over all banned states
      group_by(type, category, .draw, time) %>% 
      summarise({{target}} := sum(.data[[target]]), 
                denom = sum(.data[[denom]]), 
                ypred=sum(ypred), 
                mu = log(sum(exp(mu))),
                mu_treated = log(sum(exp(mu_treated))),
                years=mean(interval(start_date, end_date) / years(1)))

  } else {
    merged_df <- merged_df %>%
      filter(state == target_state, exposure_code == 1) %>%
      mutate(years = interval(start_date, end_date) / years(1), denom=.data[[denom]])
  }
  
  table_df <- merged_df %>%
    ungroup() %>%
    ## Aggregate over time
    group_by(type, category, .draw) %>%
    summarize(
      ypred = sum(ypred),
      outcome = sum(.data[[target]]), years = mean(years),
      treated = sum(exp(mu_treated)), untreated = sum(exp(mu)),
      denom = ifelse(target == "births", sum(denom * years, na.rm = TRUE), sum(denom, na.rm = TRUE)),
      treated_rate = treated / denom * rate_normalizer,
      untreated_rate = untreated / denom * rate_normalizer,
      outcome_rate = round(outcome / denom * rate_normalizer, 2),
      outcome_diff = round(treated - untreated)
    ) %>%
    ungroup() %>%
    ## Compute quantiles of effects
    group_by(type, category) %>%
    summarize(
      ypred_mean = mean(ypred),
      outcome = mean(outcome),
      outcome_diff_mean = round(mean(outcome_diff)), 
      outcome_diff_lower = round(quantile(outcome_diff, 0.025)), 
      outcome_diff_upper = round(quantile(outcome_diff, 0.975)),
      outcome_rate = mean(outcome_rate),
      ypred_lower = quantile(ypred, 0.025), ypred_upper = quantile(ypred, 0.975),
      treated_mean = mean(treated), treated_lower = quantile(treated, 0.025), treated_upper = quantile(treated, 0.975),
      untreated_mean = mean(untreated), untreated_lower = quantile(untreated, 0.025), untreated_upper = quantile(untreated, 0.975),      
      treated_rate_mean = mean(treated_rate), treated_rate_lower = quantile(treated_rate, 0.025), treated_rate_upper = quantile(treated_rate, 0.975),
      untreated_rate_mean = mean(untreated_rate), untreated_rate_lower = quantile(untreated_rate, 0.025), untreated_rate_upper = quantile(untreated_rate, 0.975), 
      causal_effect_diff_mean = mean(treated_rate - untreated_rate), causal_effect_diff_lower = quantile(treated_rate - untreated_rate, 0.025), causal_effect_diff_upper = quantile(treated_rate - untreated_rate, 0.975),
      causal_effect_ratio_mean = mean(treated_rate / untreated_rate), causal_effect_ratio_lower = quantile(treated_rate / untreated_rate, 0.025), causal_effect_ratio_upper = quantile(treated_rate / untreated_rate, 0.975),
      denom = mean(denom),
      pval = 2*min(mean(untreated_rate > treated_rate), mean(untreated < treated))
     )
    
  table_df <- table_df %>%
  mutate(
    # ypred_mean_rate = ypred_mean / years / (denom / rate_normalizer),
    rate_diff = round(causal_effect_diff_mean, 2),
    rate_diff_lower = round(causal_effect_diff_lower, 2),
    rate_diff_upper = round(causal_effect_diff_upper, 2),
    mult_change = causal_effect_ratio_mean,
    mult_change_lower = causal_effect_ratio_lower,
    mult_change_upper = causal_effect_ratio_upper)

  table_df <- table_df %>%
    mutate(birth_counts_str = paste0(outcome_diff_mean, " (", outcome_diff_lower, ", ", outcome_diff_upper, ")")) %>%
    mutate(birth_rate_abs_str = paste0(rate_diff, " (", rate_diff_lower, ", ", rate_diff_upper, ")")) %>%
    mutate(birth_rate_pct_str = paste0(round(100*(mult_change-1), 2), " (", round(100*(mult_change_lower-1), 2), ", ", round(100*(mult_change_upper-1), 2), ")")) %>%
    ungroup() %>%
    filter((type != "total" & category !="Total")|type == "total")
  

  pvals <- pval_rows <- table_df %>% pull(pval)
  pval_rows <- which(pvals < 0.05)
  table_df <- table_df %>% mutate(category = paste0(category, ifelse(pval <= 0.05, "*", "")))
  
  table_df %>%
    select(type, category, outcome, outcome_rate, birth_counts_str, birth_rate_abs_str, birth_rate_pct_str) %>%
    gt(rowname_col = "category") |>
  tab_header(
    title = tab_caption
  ) |> 
  ## ROW OPERATIONS
  tab_row_group(
    label = "Insurance type",
    rows = type == "insurance"
  ) |>
  tab_row_group(
    label = "Education",
    rows = type == "edu"
  ) |>
  tab_row_group(
    label = "Marital status",
    rows = type == "marital"
  ) |>
  tab_row_group(
    label = "Race/ethnicity",
    rows = type == "race"
  ) |>
  tab_row_group(
    label = "Age",
    rows = type == "age"
  ) |>
  row_group_order(groups = c(NA, "Age", "Race/ethnicity", "Marital status",
                             "Education", "Insurance type")) |>
  ### COLUMN OPERATIONS
  tab_spanner(
    label = "Birth rate (per 1,000 women per year)",
    columns = c(outcome_rate, birth_rate_abs_str, birth_rate_pct_str)) |>
  tab_spanner(
    label = "Birth count",
    columns = c(outcome, birth_counts_str)) |>
  cols_label(
    outcome_rate = "Observed",
    birth_rate_abs_str = html("Difference from expected <br>(95% CI)"),
    birth_rate_pct_str = html("Percent change from expected<br>(95% CI)"),
    outcome = "Observed",
    birth_counts_str = html("Difference from expected<br>(95% CI)"),
    category = ""
  ) |>
  tab_stub_indent(
    rows = category != "Total",
    indent = 5
  ) -> table_df
  
  ## Styling
  table_df |>
  tab_options(table.align = "left", heading.align = "left") |>
  cols_align(align = "left") |>
  cols_hide(c(type, category)) |>
  tab_options(table.font.size=8) |>
  opt_vertical_padding(scale = 0.5) |>
  cols_width(category ~ px(125),
             birth_rate_abs_str ~ px(100),
             birth_rate_pct_str ~ px(100),
             outcome_rate ~ px(50),
             birth_counts_str ~ px(100),
             outcome ~ px(50)) -> table_df_final

  table_df_final
}