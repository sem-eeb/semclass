##################
# Code for analysis of categorical data
# in SEMs
#
# Jarrett Byrnes
##################

library(piecewiseSEM)
library(nlme)

#data from https://github.com/jebyrnes/phrag_common_garden_sem
#Bowen, J.L., Kearns, P.J., Byrnes, J.E.K., Wigginton, S., Allen, W.J., 
#Greenwood, M., Tran, K., Yu, J., Cronin, J.T., Meyerson, L.A., 2017. Lineage 
#overwhelms environmental conditions in determining rhizosphere bacterial 
#community structure in a cosmopolitan invasive plant. Nature Communications 
#8, 501. 

bowen <- read.csv("../data/bowen.csv")

###
# A simple categorical model
###
div_mod <- lme(observed_otus ~  status, 
                      random =~ 1|Genotype,
                      data = bowen, method = "ML")

activity_mod <- lme(RNA.DNA ~ status + observed_otus, 
                           random =~ 1|Genotype, 
                           data=bowen, method="ML")

c_mod <- lme(below.C ~ observed_otus +  status, 
                    random =~ 1|Genotype, data=bowen, method="ML")

biomass_mod <- lme(abovebiomass_g ~ RNA.DNA + observed_otus + below.C + status, 
                           random =~ 1|Genotype, data = bowen, method="ML")

bowen_mod <- psem(
  div_mod,
  activity_mod,
  c_mod,
  biomass_mod,
  data = bowen
)

basisSet(bowen_mod)
dSep(bowen_mod)
coefs(bowen_mod)

#what do the ANOVAs say?
anova(bowen_mod)

#The coefficients
library(emmeans)

lapply(bowen_mod[-length(bowen_mod)], emmeans, specs = ~status )


###
# Test of Differences of means
###

#Let's look at, posthoc tests
generic_tukey <- function(x)
  emmeans(x, list(pairwise ~ status))


lapply(bowen_mod[-length(bowen_mod)], generic_tukey)

####
## Multigroup ####
####

meadows<-read.csv("../data/FinnishMeadows_Multigroup.csv")


###########################
# Penguin Example      ####
###########################
#data
library(palmerpenguins)

#for figures
library(ggplot2)
library(patchwork)

penguins <- na.omit(penguins)

#view some relationships
flipper_plot <- ggplot(data = penguins,
                       aes(x = body_mass_g, 
                           y = flipper_length_mm,
                           color = species, fill = species)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_minimal(base_size = 17) +
  scale_color_manual(values = c("darkorange","purple","cyan4")) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(fill = "", color = "")


size_plot <- ggplot(data = penguins,
                    aes(x = species, 
                        y = body_mass_g,
                        fill = species)) +
  geom_boxplot() +
  theme_minimal(base_size = 17) +
  scale_fill_manual(values = c("darkorange","purple","cyan4"))  +
  guides(fill = "none")


flipper_plot + size_plot + plot_layout(guides = "collect")

#fit a model with no direct species effect
size_mod <- lm(body_mass_g ~ species ,
               data = penguins)


flipper_mod <- lm(flipper_length_mm ~ body_mass_g,
                  data = penguins)

penguin_fit <- psem(
  size_mod,
  flipper_mod)

#check fit
dSep(penguin_fit)
fisherC(penguin_fit)
LLchisq(penguin_fit)

#fit a flipper model with species
flipper_mod_sp <- lm(flipper_length_mm ~ 
                       body_mass_g + species,
                     data = penguins)

penguin_fit_species <- psem(
  size_mod,
  flipper_mod_sp)

#check fit
fisherC(penguin_fit_species)
LLchisq(penguin_fit_species)
AIC(penguin_fit_species)


# interactions
flipper_mod_int <- lm(flipper_length_mm ~ 
                        body_mass_g * species,
                      data = penguins)

penguin_fit_int <- psem(
  size_mod,
  flipper_mod_int)

#compare the models
AIC(penguin_fit)
AIC(penguin_fit_species)
AIC(penguin_fit_int)

anova(penguin_fit, penguin_fit_species, penguin_fit_int)

#can we see interaction coefs?
coefs(penguin_fit_int)

#oops, no, so, emmeans!
library(emmeans)

emmeans(size_mod, ~species)

#get average mass by species
weight_list <- split(penguins$body_mass_g, penguins$species)
weight_vec <- sapply(weight_list, mean)

emmeans(flipper_mod_int, ~ species|body_mass_g, 
        at = list(body_mass_g =weight_vec) )

emtrends(flipper_mod_int, 
         ~species,
         var = "body_mass_g")

###########################
# OVB Models           ####
###########################
library(lme4)

temp_mod <- lmer(temperature ~ (1|Site), data = dat)

snail_mod <- lmer(snails ~ temperature + (1|Site), 
                  data = dat)


