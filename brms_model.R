library(brms)
library(geiger)
library(phytools)
library(ape)
library(dplyr)
library(assertthat)
library(ggplot2)
library(haven)
library(tidybayes)
library(ggridges)
library(RColorBrewer)

set.seed(1990)

if(!dir.exists('models'))
  dir.create('models')

data = readRDS("Schulzdata_withGlottocodes.rds")

# Get complete dataset
data = data[!is.na(data$lang_glottocode) & !is.na(data$ChurchExpWest) & !is.na(data$KII),]

assert_that(nrow(data) == 160, msg = "The contents of the data has changed")
jager_tree = read.tree("jaeger_pruned.tree")

## Match data to phylogeny
data = data[data$lang_glottocode %in% jager_tree$tip.label,]
jager_tree = keep.tip(jager_tree, data$lang_glottocode)

# Check data is the size we expect
assert_that(all(data$lang_glottocode %in% jager_tree$tip.label), msg = "Tree tips do not line up")

# Create a dataframe of the variables we want to model
model_df = data[,c("lang_glottocode", "ChurchExpWest", "ChurchExpEast", "KII", "Country")]
# Ensure all data have values for all variables
model_df = model_df[complete.cases(model_df),] %>% 
  as.data.frame()

# A number of glottocodes are applied to multiple nations (e.g. Spanish in SA)
# We use a repeated measures model to account for this relationship. 
# Secondly, we run a model with one nation per taxa to test the sensitivity of this analysis
model_df$species = model_df$lang_glottocode

# Tree to covariance matrix for the model
A = ape::vcv.phylo(jager_tree)

# First, we run a model replicating the relationship in Schulz et al (2019)
# i.e. without phylogenetic controls
fit_uncontrolled = brm(
  KII ~ ChurchExpWest + ChurchExpEast, 
  data = model_df, 
  family = gaussian(), 
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,  save_all_pars = TRUE,
  iter = 4000, warmup = 1000, file = "models/uncontrolled_model.rds"
)

# Model summary shows a coefficient of ChurchExpWest = -0.23
summary(fit_uncontrolled)

# Phylogenetically controlled and repeated measures model
fit_controlled <- brm(
  KII ~ ChurchExpWest + ChurchExpEast + (1|gr(lang_glottocode, cov = A)) + (1|species), 
  data = model_df, 
  family = student(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,  save_all_pars = TRUE,
  iter = 4000, warmup = 1000, file = "models/controlled_model.rds"
)

# Model summary 
summary(fit_controlled)

# Phylogenetically controlled model, but only allowing 1 society / language family
# model_dfu = model_df[!duplicated(model_df$lang_glottocode),]
model_dfu = model_df

# Dutch should only link to Netherlands - remove Belgium
model_dfu = model_dfu[!(model_dfu$Country == "Belgium" & model_dfu$lang_glottocode == "dutc1256"),]
# Fang should be Equiatorial Guinea not Gabon
model_dfu = model_dfu[!(model_dfu$Country == "Gabon" & model_dfu$lang_glottocode == "fang1246"),]
# Hungarian should be Hungary, not Austria
model_dfu = model_dfu[!(model_dfu$Country == "Austria" & model_dfu$lang_glottocode == "hung1274"),] 
# link Korean to South Korea, rather than North Korea
model_dfu = model_dfu[!(model_dfu$Country == "North Korea" & model_dfu$lang_glottocode == "kore1280"),] 
# Modern Greek linked to Greece no Cyprus
model_dfu = model_dfu[!(model_dfu$Country == "Cyprus" & model_dfu$lang_glottocode == "mode1248"),] 
# North Levantine Arabic linked to Lebannon not Syria
model_dfu = model_dfu[!(model_dfu$Country == "Syria" & model_dfu$lang_glottocode == "nort3139"),] 
# Portuguese linked to Portugal not Brazil
model_dfu = model_dfu[!(model_dfu$Country == "Brazil" & model_dfu$lang_glottocode == "port1283"),]  
# Romanian linked to Romania and not Moldova
model_dfu = model_dfu[!(model_dfu$Country == "Moldova" & model_dfu$lang_glottocode == "roma1327"),]  
# Gulf Arabic linked to United Arab Emirates and not Kuwait or Qatar. 
model_dfu = model_dfu[!(model_dfu$Country %in% c("Kuwait", "Qatar") & model_dfu$lang_glottocode == "gulf1241"),]  
# Russian is linked to the Russian Federation and not Belarus or Kazakhstan 
model_dfu = model_dfu[!(model_dfu$Country %in% c("Belarus", "Kazakhstan") & model_dfu$lang_glottocode == "russ1263"),]   
# English is linked to United Kingdom
model_dfu = model_dfu[!(model_dfu$Country %in% c("Australia", "Canada", "Ireland", "New Zealand", "United States") & model_dfu$lang_glottocode == "stan1293"),]   
# Spanish is linked to Spain
model_dfu = model_dfu[!(model_dfu$Country %in% c("Argentina", "Belize", 
                        "Chile", "Colombia", "Costa Rica", "Cuba", 
                        "Dominican Republic", "Ecuador", "Guatemala", 
                        "Honduras", "Mexico", "Nicaragua", "Panama", "Peru", 
                        "Puerto Rico", "El Salvador", "Uruguay", "Venezuela") & 
                          model_dfu$lang_glottocode == "stan1288"),]    

fit_controlledunique <- brm(
  KII ~ ChurchExpWest + ChurchExpEast + (1|gr(lang_glottocode, cov = A)), 
  data = model_dfu, 
  family = student(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0,50), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma"),
    prior(student_t(3,0,10), "nu")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,  save_all_pars = TRUE, 
  iter = 10000, warmup = 6000, file = "models/controlled_modelunique.rds",
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit_controlledunique)

### Model checks ####
### Uncontrolled model
residuals_uc = add_residual_draws(newdata = model_df, fit_uncontrolled)

# Normal QQ
residuals_uc %>%
  median_qi() %>%
  ggplot(., aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()

### Controlled model
residuals_c = add_residual_draws(newdata = model_df, fit_controlled)

# Normal QQ
residuals_c %>%
  median_qi() %>%
  ggplot(., aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()

### Unique controlled model 
residuals_cu = add_residual_draws(newdata = model_df, fit_controlledunique)

residuals_cu %>%
  median_qi() %>%
  ggplot(., aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()


# Plot predicted values against true values
predictedoutput_controlled = predict(fit_controlled)
predictedoutput_uncontrolled = predict(fit_uncontrolled)

plot(x = model_df$KII, y = predictedoutput_controlled[,1])
points(x = model_df$KII, y = predictedoutput_uncontrolled[,1], col = "red")

## Calculate pagel's lambda for controlled model
hyp <- paste(
  "sd_lang_glottocode__Intercept^2 /", 
  "(sd_lang_glottocode__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(fit_controlled, hyp, class = NULL)
hyp
plot(hyp)

# Species level effects (independent of phylogeny)
hyp <- paste(
  "sd_species__Intercept^2 /", 
  "(sd_lang_glottocode__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(fit_controlled, hyp, class = NULL)
hyp
plot(hyp)

# Calculate lambda for Unique controlled model (no repeated measures here)
hyp <- paste(
  "sd_lang_glottocode__Intercept^2 /", 
  "(sd_lang_glottocode__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(fit_controlledunique, hyp, class = NULL)
hyp
plot(hyp)

### Model comparisons
fit_uncontrolled = add_criterion(fit_uncontrolled, "loo", moment_match = TRUE)
fit_controlled = add_criterion(fit_controlled, "loo", moment_match = TRUE)
fit_controlledunique = add_criterion(fit_controlledunique, "loo", moment_match = TRUE)

bayes_R2(fit_uncontrolled)
bayes_R2(fit_controlled)

# Can only comapre uncontrolled and controlled, because unique has a different dataset
loo_compare(fit_uncontrolled, fit_controlled)

### Plot estimates of coefficent
post_uncontrolled = posterior_samples(fit_uncontrolled)
post_controlled = posterior_samples(fit_controlled)

plot_df = data.frame(value = c(post_controlled$b_ChurchExpWest, 
                               post_uncontrolled$b_ChurchExpWest),
                     model = rep(c("Phylogenetically controlled", "No phylogenetic control"), each = 6000))


pdf("model_comparison.pdf", width = 7, height = 5)
ggplot(plot_df, aes(x = value, y = model, fill = factor(stat(quantile)))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = c(0.025, 0.975), quantile_lines = TRUE
  )  + 
  scale_fill_manual(
    name = "Probability", values = c("#A0A0A0A0", "#FF0000A0", "#A0A0A0A0"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  ) + 
  geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("") + xlab("") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  ggtitle("Beta Coefficient for Western Church exposure", "Posterior distribution with 95% intervals") 
dev.off()
  

#### Phylo Signal ####
kii = model_dfu$KII
names(kii) = model_dfu$lang_glottocode
church = model_dfu$ChurchExpWest
names(church) = model_dfu$lang_glottocode

phylosig(jager_tree, kii, method = "lambda")
phylosig(jager_tree, church, method = "lambda")
