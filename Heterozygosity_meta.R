## JClarke alt data
library(brms)
library(lme4)
library(tidyverse)
library(readxl)
library(cmdstanr)



log_ratio <- function(x1,x2){
  rt <- log(x2/x1)
  return(rt)
}

# this approximation of the ratio of two normal variates =
# equation 8 in  https://doi.org/10.1007/s00362-012-0429-2 
se_log_ratio <- function(x1, #ratio denominator
                         x2, #ratio numerator
                         se1, #se of denominator
                         se2){ # se of numerator){
  
ro = (se2/se1)^-2
beta = (x2/x1)^2

se_rat <- sqrt(se2^2*(ro + beta))
rt_se_rat <- log((x2/x1) + se_rat)-log(x2/x1) # approximate se of log ratio as difference between log(ratio + 1SE)-log(ratio)

  return(rt_se_rat)
  
}


# Julia's code ------------------------------------------------------------

p_positive <- function(x){
  length(which(x > 0))/length(x)
}
df_alt <- read_xlsx("He_dGR.xlsx", sheet = 1) %>%
  select(Study:fit_par) %>% 
  mutate(HR = Post_He/Pre_He,
         dHe = log_ratio(Pre_He,Post_He),
         dHE_er = se_log_ratio(Pre_He,Post_He,PreHe_SE,PostHe_SE),
         GR = Post_fit/Pre_fit,
         dGR = log_ratio(Pre_fit,Post_fit),
         smd_GR_er = se_log_ratio(Pre_fit,Post_fit,Prefit_SE,Postfit_SE),
         Study = factor(Study),
         ID = row_number())  %>% 
  drop_na()

mean(df_alt$HR)

tst <- ggplot(data = df_alt,
              aes(x = dHe,y = dGR))+
geom_point()+
  # geom_errorbar(aes(ymin = dGR-smd_GR_er,
  #                   ymax = dGR+smd_GR_er))+
  geom_errorbarh(aes(xmin = dHe-dHE_er,
                     xmax = dHe+dHE_er))+
  geom_smooth(method = "lm")
  
tst

hist(df_alt$PreHe_SE,
     breaks = 30)


tst <- ggplot(data = df_alt,
              aes(x = Pre_He,y = PreHe_SE))+
  geom_point()+
  geom_smooth(method = "lm")

tst



# custom Stan model -------------------------------------------------------

n_studies <- as.integer(max(df_alt$Study_n))

stan_data <- list(
  n_obs = nrow(df_alt), #number of observations
  He_pre = df_alt$Pre_He, # vector of pre He values
  He_post = df_alt$Post_He, #vector of post He values
  He_se_pre = df_alt$PreHe_SE, #vector of pre SE of He values
  He_se_post = df_alt$PostHe_SE, #vector of post SE of He values
  
  n_studies = n_studies, # number of studies 
  study = df_alt$Study_n # vector of study indicators
  
  
)


mod_f <- "models/stan_measure_error.stan"

mod <- cmdstanr::cmdstan_model(mod_f)

fit <- mod$sample(
  data = stan_data,
  parallel_chains = 4,
  iter_sampling = 5000,
  iter_warmup = 5000,
  adapt_delta = 0.95
)

summ <- posterior::as_draws_df(fit) %>% 
  posterior::summarise_draws("mean",
                             "sd",
                             "ess_bulk",
                             "rhat",
                             ~ quantile(.x,probs = c(0.025,0.05,0.1,0.25,0.75,0.9,0.95,0.975)),
                             p_positive = p_positive) %>% 
  rename_with(.,.cols = contains("%"),
              ~ paste0("CI_",gsub(".","_",gsub("%","",.x), fixed = TRUE))) 

# summ <- fit$summary()

beta_fit <- summ %>% 
  filter(grepl("beta[",variable, fixed = TRUE)) 

BETA_fit <- summ %>% 
  filter(grepl("BETA",variable, fixed = TRUE)) %>% 
  mutate(Study = "Overall mean")

beta_out <- df_alt %>% 
  select(Study,Study_n) %>% 
  arrange(Study_n) %>%
  distinct() %>% 
  bind_cols(.,beta_fit) %>% 
  arrange(mean) %>% 
  bind_rows(.,BETA_fit)

write_csv(beta_out,file = "output/differences_heterozygosity.csv")

beta_out <- beta_out %>% 
  mutate(Study_sort = factor(Study,
                             levels = beta_out$Study,
                             ordered = TRUE))

BETA_plot <- beta_out %>% 
  filter(Study == "Overall mean")
delta_he <- ggplot(data = beta_out)+
  geom_errorbar(aes(x = Study_sort,
                    y = mean,
                    ymin = CI_2_5,
                    ymax = CI_97_5),
                width = 0,
                alpha = 0.3)+
  geom_errorbar(aes(x = Study_sort,
                    y = mean,
                    ymin = CI_10,
                    ymax = CI_90),
                width = 0,
                alpha = 0.5,
                linewidth = 1)+
  geom_errorbar(data = BETA_plot,
                aes(x = Study_sort,
                    y = mean,
                    ymin = CI_2_5,
                    ymax = CI_97_5),
                width = 0,
                alpha = 0.8)+
  geom_errorbar(data = BETA_plot,
                aes(x = Study_sort,
                    y = mean,
                    ymin = CI_10,
                    ymax = CI_90),
                width = 0,
                alpha = 1,
                linewidth = 1.1)+
  geom_point(aes(x = Study_sort,
                 y = mean))+
  geom_hline(yintercept = 0)+
  coord_flip()+
  xlab("")+
  ylab("Predicted change in heterozygosity")+
  theme_bw()+
  theme(text = element_text(family = "serif"))


pdf("Figure_1.pdf",
    width = 6.5,
    height = 5)
print(delta_he)
dev.off()







me_prior <- c(prior(normal(0,1),
                      class = "meanme"),
              prior(normal(0,1),
                    class = "sdme"),
              prior(normal(0,1),
                    class = "b"))


#|se(smd_GR_er, sigma = TRUE)
mod4 <- brm(dGR ~ me(dHe,dHE_er) + (1|Study),     #|se(smd_GR_er, sigma = TRUE)
            family = "gaussian", 
            data = df_alt,
            prior = me_prior,
            iter = 15000,
            warmup = 10000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.95),
            backend = "cmdstanr")



summary(mod4)
pp_check(mod4,ndraws = 100)

#|se(smd_GR_er, sigma = TRUE)
mod5 <- brm(dGR ~ me(dHe,dHE_er) + tran_R + (1|Study), 
            family = "gaussian", 
            data = df_alt,
            prior = me_prior,
            iter = 10000,
            warmup = 8000,
            control = list(adapt_delta = 0.95),
            backend = "cmdstanr")

summary(mod5)

#|se(smd_GR_er, sigma = TRUE)
mod6 <- brm(dGR ~ me(dHe,dHE_er) + tran_R + fit_par + (1|Study), 
            family = "gaussian", 
            data = df_alt,
            prior = me_prior,
            iter = 10000,
            warmup = 8000,
            control = list(adapt_delta = 0.95),
            backend = "cmdstanr")

summary(mod6)


draws_4 <- posterior::as_draws_df(mod4) %>% 
  posterior::summarise_draws("mean",
                             "sd",
                             "ess_bulk",
                             "rhat",
                             ~ quantile(.x,probs = c(0.025,0.05,0.1,0.25,0.75,0.9,0.95,0.975)),
                             p_positive = p_positive) %>% 
  rename_with(.,.cols = contains("%"),
              ~ paste0("CI_",gsub(".","_",gsub("%","",.x), fixed = TRUE)))  %>% 
  mutate(model = "dHe_27")

draws_5 <- posterior::as_draws_df(mod5) %>% 
  posterior::summarise_draws("mean",
                             "sd",
                             "ess_bulk",
                             "rhat",
                             ~ quantile(.x,probs = c(0.025,0.05,0.1,0.25,0.75,0.9,0.95,0.975)),
                             p_positive = p_positive) %>% 
  rename_with(.,.cols = contains("%"),
              ~ paste0("CI_",gsub(".","_",gsub("%","",.x), fixed = TRUE)))  %>% 
  mutate(model = "dHe_tranR_27")

draws_6 <- posterior::as_draws_df(mod6) %>% 
  posterior::summarise_draws("mean",
                             "sd",
                             "ess_bulk",
                             "rhat",
                             ~ quantile(.x,probs = c(0.025,0.05,0.1,0.25,0.75,0.9,0.95,0.975)),
                             p_positive = p_positive) %>% 
  rename_with(.,.cols = contains("%"),
              ~ paste0("CI_",gsub(".","_",gsub("%","",.x), fixed = TRUE))) %>% 
  mutate(model = "dHe_tranR_fitPar_27")


draws_summary <- bind_rows(draws_4,
                           draws_5,
                           draws_6)
write_csv(draws_summary,
          "output/draws_summary.csv")

pred_data <- data.frame(dHe = seq(min(df_alt$dHe),max(df_alt$dHe),
                                  length.out = 100),
                        dHE_er = median(df_alt$dHE_er),
                        smd_GR_er = median(df_alt$smd_GR_er))
preds_4 <- predict(mod4,
                    newdata = pred_data,
                    re_formula = NA,
                    probs = c(0.025,0.075,0.925,0.975))
preds_4 <- bind_cols(pred_data,
                      preds_4)
preds_4_line <- preds_4[c(1,nrow(preds_4)),]


comp_plot <- ggplot()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_ribbon(data = preds_4,
              aes(x = dHe, ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.1)+
  geom_ribbon(data = preds_4,
              aes(x = dHe, ymin = Q7.5,
                  ymax = Q92.5),
              alpha = 0.1)+
  geom_line(data = preds_4_line,
            aes(x = dHe, y = Estimate))+
  geom_errorbarh(data = df_alt,
                 aes(y = dGR,xmin = dHe-dHE_er,
                     xmax = dHe+dHE_er),
                 alpha = 0.2)+
  geom_point(data = df_alt,
             aes(x = dHe,y = dGR, colour = fit_par,
                 size = tran_R))+
  # scale_colour_brewer(type = "qual",
  #                    palette = "Set3")+
  scale_colour_viridis_d(begin = 0,end = 0.9,
                         direction = -1)+
  theme_bw()+
  ylab("log(GR)")+
  xlab("log(HR)")+
  theme(legend.position = "right")
print(comp_plot)

pdf("Figure_2.pdf",
    width = 6.5,
    height = 6)
print(comp_plot)
dev.off()

refit_loo <- FALSE # change to true to re-run cross-validation
if(refit_loo){ # within this conditional block because requires many hours to run
  
  exact_loo6 <- kfold(mod6,folds = "loo")
  exact_loo4 <- kfold(mod4,folds = "loo")
  exact_loo5 <- kfold(mod5,folds = "loo")
  

loo_out <- loo_compare(exact_loo4,
                        exact_loo5,
                        exact_loo6)
loo_out <- as_tibble(loo_out) %>% 
  mutate(model = rownames(loo_out))


write_csv(loo_out,
          "output/loo_out.csv")
}

modls <- data.frame(model = c("mod4","mod5","mod6"),
                    model2 = c("dHe_27","dHe_tranR_27","dHe_tranR_fitPar_27"))
param_HR <- draws_summary %>% 
  filter(variable == "bsp_medHedHE_er") %>% 
  select(mean,CI_2_5,CI_97_5,p_positive,model) %>% 
  mutate(p_positive = p_positive*100,
         across(.cols = mean:p_positive, .fns = ~signif(.x,2)),
         HR_effect = paste0(mean," [",CI_2_5," : ",CI_97_5,"] ",p_positive,"% positive")) %>% 
  select(model,HR_effect)

param_tran <- draws_summary %>% 
  filter(variable == "b_tran_R") %>% 
  select(mean,CI_2_5,CI_97_5,p_positive,model) %>% 
  mutate(p_positive = p_positive*100,
         across(.cols = mean:p_positive, .fns = ~signif(.x,2)),
         TR_effect = paste0(mean," [",CI_2_5," : ",CI_97_5,"] ",p_positive,"% positive")) %>% 
  select(model,TR_effect)

params <- param_HR %>% 
  left_join(.,param_tran,
            by = "model")
table_3 <- loo_out %>% 
  select(elpd_diff,se_diff,model) %>% 
  left_join(.,modls,
            by = "model") %>% 
  left_join(.,params,
            by = c("model2" = "model")) %>% 
  select(-model)

write_csv(table_3,"output/Table_3.csv")
save.image("temp.RData")
