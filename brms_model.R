library(brms)
library(geiger)
library(ape)
library(dplyr)
library(assertthat)
library(ggplot2)
library(haven)
library(tidybayes)
library(ggridges)
library(RColorBrewer)

set.seed(1990)
dir.create("models")

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
model_df = data[,c("lang_glottocode", "ChurchExpWest", "KII")]
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
  KII ~ ChurchExpWest, 
  data = model_df, 
  family = gaussian(), 
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2, criterion = "loo",
  iter = 4000, warmup = 1000, file = "models/uncontrolled_model.rds"
)

# Model summary shows a coefficient of ChurchExpWest = -0.21
summary(fit_uncontrolled)

# Phylogenetically controlled and repeated measures model
fit_controlled <- brm(
  KII ~ ChurchExpWest + (1|gr(lang_glottocode, cov = A)) + (1|species), 
  data = model_df, 
  family = student(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2, criterion = "loo",
  iter = 4000, warmup = 1000, file = "models/controlled_model.rds"
)

# Model summary 
summary(fit_controlled)

# Phylogenetically controlled model, but only allowing 1 society / language family
model_dfu = model_df[!duplicated(model_df$lang_glottocode),]

fit_controlledunique <- brm(
  KII ~ ChurchExpWest + (1|gr(lang_glottocode, cov = A)), 
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
  sample_prior = TRUE, chains = 2, cores = 2, criterion = "loo",
  iter = 10000, warmup = 6000, file = "models/controlled_modelunique.rds",
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit_controlledunique)
### model checks

residuals = add_residual_draws(newdata = model_df, fit_controlled)

ggplot(residuals, aes(x = .row, y = .residual)) +
  stat_pointinterval()

residuals %>%
  median_qi() %>%
  ggplot(., aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()


predicted_output = predict(fit_controlled)

cor(model_df$KII, predicted_output[,1])

## Calculate pagel's lambda
hyp <- paste(
  "sd_lang_glottocode__Intercept^2 /", 
  "(sd_lang_glottocode__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(fit_controlled, hyp, class = NULL)
hyp
plot(hyp)

hyp <- paste(
  "sd_species__Intercept^2 /", 
  "(sd_lang_glottocode__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(fit_controlled, hyp, class = NULL)
hyp
plot(hyp)

### Model comparisons

bayes_R2(fit_uncontrolled)
bayes_R2(fit_controlled)

loo_compare(fit_uncontrolled, fit_controlled)

### Plot estimates of coefficent

post_uncontrolled = posterior_samples(fit_uncontrolled)
post_controlled = posterior_samples(fit_controlled)

plot_df = data.frame(value = c(post_controlled$b_ChurchExpWest, 
                               post_uncontrolled$b_ChurchExpWest),
                     model = rep(c("Phylogenetically controlled", "No phylogenetic control"), each = 6000))


pdf("model_comparison.pdf")
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
    
### Plot tree

model_df$colours = brewer.pal(9, "YlOrRd")[cut(model_df$ChurchExpWest, breaks = 8)]

tree = ladderize(jager_tree)

plot(tree, show.tip.label = FALSE)
tiplabels(pch=21, bg=(model_df$colours),cex=1, )
