library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(janitor)
library(surveillance)
library(broom)
library(jsonlite)

theme_set(theme_light())
## Setup a ggplot theme to be used in the report.
our_facet_theme <- theme(strip.background = element_rect(fill = "#155540"),
        strip.text = element_text(colour = 'white', face = "bold"),
        legend.position = "inside")
        
#file_name <- here::here("data", "SN.csv")
#df <- vroom::vroom(file_name) %>% 
# mutate_at(vars(contains("date")), dmy)

now <- as.POSIXlt(Sys.time())
last_valid_date <- now - days(2) # Meldeverzug

data <- fromJSON(paste0('https://services7.arcgis.com/mOBPykOjAyBO2ZKk/arcgis/rest/services/RKI_COVID19/FeatureServer/0/query?where=%28neuerfall%3D0+or+neuerfall%3D1%29+and+idbundesland%3D14+and+refdatum%3C%3Emeldedatum+and+meldedatum%3C\'', format(last_valid_date, "%Y-%m-%d"), '\'&objectIds=&time=&resultType=none&outFields=Refdatum%2CMeldedatum&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnDistinctValues=false&cacheHint=false&orderByFields=Refdatum&groupByFieldsForStatistics=Refdatum&outStatistics=&having=&resultOffset=&resultRecordCount=&sqlFormat=none&f=pjson&token='))

df2 <- flatten(data$features)
df2 <- transmute(df2, date_confirmation = format(as.Date(as.POSIXct(`attributes.Meldedatum` / 1000, origin="1970-01-01")), "%d.%m.%Y"), date_onset_symptoms = format(as.Date(as.POSIXct(`attributes.Refdatum` / 1000, origin="1970-01-01")), "%d.%m.%Y"))
df <- df2 %>% mutate_at(vars(contains("date")), dmy)

delay_cutoff <- 21 # days


impute_delay <- function() {
  
  train <- df %>% 
    mutate(delay = as.numeric(date_confirmation - date_onset_symptoms)) %>% 
    filter(delay >= 0 ) %>% 
    mutate(week_confirm = week(date_confirmation)) %>% 
    mutate(delay_prep_log = if_else(delay == 0, 1e-2, delay),
           event_observed = 1) 
  
  #Compute summary of the delay distribution
  delay_summary <- train %>% 
    group_by(week_confirm) %>%
    summarise(n = n(), 
              n_delay_available = sum(!is.na(delay)),
              prop_delay_available = scales::percent(mean(!is.na(delay))),
              mean = mean(delay, na.rm = TRUE),
              median = median(delay, na.rm = TRUE),
              q25 = quantile(delay, prob=0.25, na.rm = TRUE),
              q75 = quantile(delay, prob=0.75, na.rm = TRUE))
  
  # Make a dataset without NAs
  train_complete <- train %>% 
    select(delay_prep_log, week_confirm, date_confirmation) %>% 
    na.omit()
  
  
  # Fit model - for some reason gamlss has problems if the dataset is not declared as a global variable?
  train_complete <<- train_complete
  impute_model <- gamlss::gamlss(
    delay_prep_log  ~ gamlss::cs(week_confirm), #as.factor(week_confirm), 
    sigma.formula = ~ gamlss::cs(week_confirm), #as.factor(week_confirm),
    family = gamlss.dist::WEI2, 
    data = train_complete, 
    method = mixed(80,80)
  )
  
  # Impute delay. For some reason this df should not contain irrelevant columns
  train_simple <- train %>% 
    select(week_confirm)
  
  # Fix seed for imputation.
  set.seed(1234)
  train_pred <- train %>% 
    mutate(
      mu = exp(predict(impute_model, newdata = train_simple, what = "mu")),
      sigma = exp(predict(impute_model, newdata = train_simple, what = "sigma")),
      delay_imputed = gamlss.dist::rWEI2(nrow(train_simple), mu = mu, sigma = sigma)
    ) 
  
  # Show model fit and empirical distirbution for each week_confirm. 
  # Put this into a function instead!
  step_width <- 0.01
  max_delay <- max(train$delay, na.rm=TRUE)
  pdf_delay <- tibble(delay = seq(0, max_delay, by = step_width)) %>% 
    mutate(delaym1 = lag(delay)) %>% 
    slice(-1) 
  
  #Range of the weeks
  week_range <- week(min(train %>% pull(date_confirmation)) ):(week((max(train %>% pull(date_confirmation)))))
  
  pdf_delay <- pdf_delay %>% crossing(data.frame(week_confirm = week_range)) %>%
    as.data.frame()

  pdf_delay <- pdf_delay %>% 
    mutate(
      mu = exp(predict(impute_model, newdata= pdf_delay %>% select(week_confirm), what="mu")),
      sigma = exp(predict(impute_model, newdata= pdf_delay %>% select(week_confirm), what="sigma")),
      prob = gamlss.dist::pWEI2(delay, mu=mu, sigma=sigma) - 
        gamlss.dist::pWEI2(delaym1, mu=mu, sigma=sigma),
      proportion_reported = gamlss.dist::pWEI2(delay, mu=mu, sigma=sigma),
      cdf = gamlss.dist::pWEI2(delay, mu=mu, sigma=sigma))
  
  ##Extract median delay from model
  model_delays <- pdf_delay %>% group_by(week_confirm) %>% do({
    q_vec <- c(0.25,0.5,0.75)
    idx <- map_int(q_vec, function(q) {
      which.min(abs(.$proportion_reported - q))
    })
    data.frame(week_confirm=.$week_confirm[idx], model_q=.$delay[idx], q=str_c("model_",q_vec))
  }) %>%  spread(q, model_q)
  
  delay_summary <- inner_join(delay_summary, model_delays, by="week_confirm")
  
  ## Show match using the cumulative distribution 
  g_cdf <- train %>% 
    arrange(week_confirm) %>% 
    mutate(week_confirm_chr = glue::glue("calendar\nweek {week_confirm}"),
           week_confirm_chr = factor(week_confirm_chr, levels = unique(week_confirm_chr))) %>% 
    ggplot(aes(x = delay)) + 
    facet_grid(. ~ week_confirm_chr) + 
    stat_ecdf() +
    ylab("Proportion reported") + xlab("Delay (days)") +
    scale_y_continuous(labels=scales::percent) + 
    geom_line(
      data = pdf_delay %>% 
        mutate(week_confirm_chr = glue::glue("calendar\nweek {week_confirm}"),
               week_confirm_chr = factor(week_confirm_chr, levels = unique(week_confirm_chr))), 
      aes(x = delay, y = cdf), size = 1.2, color = "#ff7f00") +
    
    scale_fill_manual(values = "#377eb8", guide = "none") +
    scale_x_continuous(breaks = c(0, 10, 20)) +
    our_facet_theme
    
  ##Show pdf of fit
  g_pdf <- train %>% 
    arrange(week_confirm) %>% 
    mutate(week_confirm_chr = glue::glue("calendar\nweek {week_confirm}"),
           week_confirm_chr = factor(week_confirm_chr, levels = unique(week_confirm_chr))) %>% 
    ggplot(aes(x = delay, y = ..density.., fill = "light")) + 
    geom_histogram(breaks = seq(0, max_delay, by = 1)) + 
    xlab("Delay between onset symptoms and confirmation (days)") + 
    facet_grid(. ~ week_confirm_chr) +
    ylab("Density") +
    geom_line(
      data = pdf_delay %>% 
        mutate(week_confirm_chr = glue::glue("calendar\nweek {week_confirm}"),
               week_confirm_chr = factor(week_confirm_chr, levels = unique(week_confirm_chr))), 
      aes(x = delay, y = prob * 1/step_width), size = 1.2, color = "#ff7f00") + 
    
    scale_fill_manual(values = "#377eb8", guide = "none") +
    scale_x_continuous(breaks = c(0, 10, 20)) +
    our_facet_theme

  #Attach imputation model as an attribute to the completed dataset
  attr(train_pred, "model") <- impute_model
  
  #Done, return imputed data and two graphs
  tab_g1 <- list(imputed_linelist=train_pred, plot_pdf=g_pdf, plot_cdf=g_cdf, delay_summary=delay_summary)
  return(tab_g1)
}

do_nowcast <- function(country_name, imputed_weibull) {

  #Select linelist for the country
  country_weibull <- imputed_weibull %>%
    as.data.frame()
  
  #Control variables of the nowcast - only do nowcasts for the last max_delay days
  now <- max(country_weibull$date_confirmation) 
  max_delay <- 21
  safePredictLag <- 0
  so_range <- c(min(country_weibull$date_onset_symptoms, na.rm = TRUE), now)

  #Fix nowcast time points so they don't depend on the imputed data.
  nowcastDates <- seq(from = now - safePredictLag - max_delay, to = now - safePredictLag, by = 1)
 
  ##Vector of possible change points (3) to be used in the bayes.trunc.ddcp model
  cp_vec <- seq(now, length.out=3, by="1 week")
  
  sts <- linelist2sts(
    as.data.frame(country_weibull), 
    dateCol="date_onset_symptoms", 
    aggregate.by = "1 day", 
    dRange = so_range)

  nc.control <- list(
    N.tInf.max = 4e3,
    N.tInf.prior = structure("poisgamma",
                             mean.lambda = mean(observed(sts)),
                             var.lambda = 5*var(observed(sts))
                             ),
    ##compute predictive distribution as wel, which is needed for some of the
    ##animations.
    predPMF = TRUE,
    dRange = so_range)

    nc <- nowcast(now = now, when = nowcastDates, data = as.data.frame(country_weibull),
                  dEventCol = "date_onset_symptoms",
                  dReportCol = "date_confirmation",
                  aggregate.by = "1 day",
                  D = delay_cutoff, # adjust cases up to 2 weeks back.
                  # # Assume constant delay distribution, but only within the last m=14 days
                  method = "bayes.trunc",
                  m = delay_cutoff, #only use last 14 days for the delay estimation
                  control = nc.control
                  # # Use the discrete time survival model with change-points
                  # method = "bayes.trunc.ddcp",
                  # control = modifyList(nc.control, list(
                  #   ddcp = list(
                  #     ddChangepoint = cp_vec,
                  #     logLambda = "iidLogGa",
                  #     eta.mu=rep(0,length(cp_vec)),
                  #     eta.prec=diag(rep(1,length(cp_vec))))
                  # ))
                  
    )
    
    ##Convert to tibble (in wide format)
    nc_tidy <- nc %>% 
      tidy() %>% 
      as_tibble() %>% 
      # Add prediction interval
      mutate(pi_lower = nc@pi[,,1],  pi_upper=  nc@pi[,,2]) %>% 
      # Return only time points which were nowcasted.
      filter(epoch %in% nowcastDates) %>% 
      # Restrict to relevant columns
      select(date, observed, predicted = upperbound, predicted_lower= pi_lower, predicted_upper = pi_upper )

    # Attach nowcast object 
    attr(nc_tidy, "original_stsNC") <- nc
    
    return(nc_tidy)
}

impute_delay()[[3]]

ncs <- map(
  "Japan", 
  do_nowcast, 
  imputed_weibull = impute_delay()[[1]] %>% 
    mutate(country = "Japan")
  ) %>% setNames("Japan")
ncs_df <- ncs %>% bind_rows() 

# ncs <- map( 
#   impute_delay()[[1]]  %>% mutate(country = "Japan"),
#   do_nowcast
#       )
# ncs_df <- ncs %>% bind_rows() 

# Reduce nowcast objects to only showing results during last 2 weeks
# A consequence of using D=14 is that older delays than 14 days are not 
# adjusted at all.
ncs_clean <- ncs_df %>% 
  filter(date > (max(date) - days(delay_cutoff))) %>% 
  mutate(obnyr = predicted - observed) %>% 
  select(date, observed, obnyr) %>% 
  gather(key, value, -date) %>% 
  ungroup()

# last case of linelist to filter nowcast data
last_linelist_case <- df %>% 
  summarise(n = n(),
            last_case = max(date_confirmation)) %>% 
  pull(last_case)

png(filename="faithful.png", width=1280, height=960)

# Plot nowcasts with corresponding prediction intervals
ncs_clean %>% 
  mutate(key = case_when(key == "obnyr" ~ "nowcast",
                         TRUE ~ key)) %>% 
  filter(date <= ymd(last_linelist_case) - days(1)) %>% 
 
  ggplot(aes(date, value)) +
  geom_col(aes(fill = key), alpha = 0.9)  +
  geom_errorbar(
    data = ncs_df %>% 
      mutate(value=1, key=NA) %>%
      filter(date > (max(date) - days(delay_cutoff)) & date < (max(date) - days(0))), 
    aes(ymin = predicted_lower, ymax = predicted_upper), width = 0.2, size = 1) +
#  facet_wrap(~country, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("#eb3a21", "#3e5ddf")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_x_date(date_breaks = "4 days", date_labels = "%b %d") +
  labs(x = "Erkankungsdatum",
       y = "Tägliche Neuerkrankungen mit bekanntem Erkrankungsdatum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  our_facet_theme

dev.off()

