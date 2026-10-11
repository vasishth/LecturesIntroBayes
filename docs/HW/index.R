## ----setup, include=FALSE-------------------------------------------------------------

knitr::opts_chunk$set(cache = TRUE, autodep = TRUE,
                      fold = FALSE,
                      class.source = "fold-show")

# same cutoff as print_warning
cut_warning <- function(wrn) {
  strwrap( paste0(wrn, collapse=""), 60, simplify = FALSE)[[1]]
}


## ----functions, include=FALSE, echo =FALSE--------------------------------------------
# USEFUL FUNCTIONS:
# makes a separated by commas list from a vector
list_str <- function(vector, and = TRUE, s = "`"){
    vector <- paste0(s,vector,s)
    if(and) {
        paste0(vector,collapse = ", ")
    } else {
        paste0(paste0(vector[-length(vector)],collapse = ", "), ", and ", vector[length(vector)])
    }
}
list_code <- function(vector){
    toString(paste0('"', vector,'"'))
}


## ----load, cache = FALSE, message = FALSE, warning=FALSE------------------------------
library(MASS)
## be careful to load dplyr after MASS
library(dplyr)
library(tidyr)
library(purrr)
library(extraDistr)
library(ggplot2)
library(loo)
library(bridgesampling)
library(brms)
library(bayesplot)
library(tictoc)
library(hypr)
library(bcogsci) ## need to install this via github.
library(lme4)
library(rstan)
library(SHELF)
library(rootSolve)

## Save compiled models:
rstan_options(auto_write = FALSE)
## Parallelize the chains using all the cores:
options(mc.cores = parallel::detectCores())
# To solve some conflicts between packages:
select <- dplyr::select
extract <- rstan::extract


## install using devtools::install_version("SIN", version= "0.6.0")
library(SIN)
#this can installed regularly
library(papaja)
library(grid)
library(gridExtra)

library(kableExtra)
library(cowplot)
library(pdftools)

## Force to load SHELF, so that the ggplto theme is active
## suppressWarnings(dist_cons <- fitdist(vals=0,probs=0))
## plotfit(dist_cons, returnPlot = TRUE)


rstan_options(silent = TRUE,
              open_progress = FALSE,
              show_messages = FALSE)
## Defauls values of some parameters
# I don't want windows opening:
formals(stan)$open_progress <- FALSE
formals(stan)$chains <- 4

## Look and feel:
# Plots
set_plots<-  function(){
bayesplot_theme_set(theme_light())
theme_set(theme_light())

# center all the plot titles?
theme_update(plot.title = element_text(hjust = 0.5))

if(knitr::is_html_output()) {
## Good for black and white and color
  options(ggplot2.continuous.colour = scale_color_viridis_c)
  options(ggplot2.continuous.fill = scale_fill_viridis_c)
  options(ggplot2.discrete.colour = scale_color_viridis_d)
  options(ggplot2.discrete.fill = scale_fill_viridis_d)
  color_scheme_set("viridis")
}

if(knitr::is_latex_output()){
  scale_color_grey_c <- function(...)  scale_colour_gradient(..., low = "white", high = "black")
  scale_fill_grey_c <- function(...)  scale_fill_gradient(..., low = "white", high = "black")
  options(ggplot2.discrete.colour = scale_color_grey)
  options(ggplot2.discrete.fill = scale_fill_grey)
  options(ggplot2.continuous.colour =  scale_color_grey_c)
  options(ggplot2.continuous.fill = scale_fill_grey_c)
  color_scheme_set("gray")
}
}
set_plots()
## IMPORTANT: this needs to be set up again in priors chapter after the functions from SHELF are used



## options(ggplot2.discrete.fill = function(...) scale_fill_viridis_d(..., option = "inferno"))

# format
options(
    htmltools.dir.version = FALSE,
    formatR.indent = 2,
    width = 70,
    digits = 3,
    signif = 2,
    warnPartialMatchAttr = FALSE,
    warnPartialMatchDollar = FALSE,
    # Don't use scientific notation:
    scipen=10,
    # tibbles:
    tibble.width = Inf,
    tibble.print_max = 3,
    tibble.print_min = 3,
    dplyr.summarise.inform = FALSE,
  tinytex.clean = FALSE #keeps the aux file for xr package
)
# output:
fixef.brmsfit <- function(...) {
  brms:::fixef.brmsfit(...) %>%
    round(2)
}
# print warning,
# last warning should be saved, e.g., saveRDS(names(last.warning),"dataR/fit_mix_rt_w.RDS")
print_warning <- function(wrn) {
cat(paste0(map_chr(strwrap(paste("Warning:", wrn), 70, simplify = FALSE), ~ paste0(.x, collapse ="\n")),"\n\n"))
}
short_summary <- function (x, digits = 2, ...)
{
  x<- summary(x)
  cat("...\n")
    # cat(" Family: ")
    # cat(summarise_families(x$formula), "\n")
    # cat("  Links: ")
    # cat(summarise_links(x$formula, wsp = 9), "\n")
    # cat("Formula: ")
    # print(x$formula, wsp = 9)
    # cat(paste0("   Data: ", x$data_name, " (Number of observations: ",
        # x$nobs, ") \n"))
    if (!isTRUE(nzchar(x$sampler))) {
        cat("\nThe model does not contain posterior samples.\n")
    }
    else {
        final_samples <- ceiling((x$iter - x$warmup)/x$thin *
            x$chains)
        # cat(paste0("Samples: ", x$chains, " chains, each with iter = ",
        #     x$iter, "; warmup = ", x$warmup, "; thin = ", x$thin,
        #     ";\n", "         total post-warmup samples = ", final_samples,
        #     "\n\n"))
        if (nrow(x$prior)) {
            cat("Priors: \n")
            print(x$prior, show_df = FALSE)
            cat("\n")
        }
        if (length(x$splines)) {
            cat("Smooth Terms: \n")
            brms:::print_format(x$splines, digits)
            cat("\n")
        }
        if (length(x$gp)) {
            cat("Gaussian Process Terms: \n")
            brms:::print_format(x$gp, digits)
            cat("\n")
        }
        if (nrow(x$cor_pars)) {
            cat("Correlation Structures:\n")
            brms:::print_format(x$cor_pars, digits)
            cat("\n")
        }
        if (length(x$random)) {
            cat("Group-Level Effects: \n")
            for (i in seq_along(x$random)) {
                g <- names(x$random)[i]
                cat(paste0("~", g, " (Number of levels: ", x$ngrps[[g]],
                  ") \n"))
                brms:::print_format(x$random[[g]], digits)
                cat("\n")
            }
        }
        if (nrow(x$fixed)) {
            cat("Population-Level Effects: \n")
            brms:::print_format(x$fixed, digits)
            cat("\n")
        }
        if (length(x$mo)) {
            cat("Simplex Parameters: \n")
            brms:::print_format(x$mo, digits)
            cat("\n")
        }
        if (nrow(x$spec_pars)) {
            cat("Family Specific Parameters: \n")
            brms:::print_format(x$spec_pars, digits)
            cat("\n")
        }
        if (length(x$rescor_pars)) {
            cat("Residual Correlations: \n")
            brms:::print_format(x$rescor, digits)
            cat("\n")
        }
        # cat(paste0("Samples were drawn using ", x$sampler, ". "))
        if (x$algorithm == "sampling") {
            #cat(paste0("For each parameter, Bulk_ESS\n", "and Tail_ESS are effective sample size measures, ",
             #   "and Rhat is the potential\n", "scale reduction factor on split chains ",
              #  "(at convergence, Rhat = 1)."))
        }
        cat("...\n")
    }
    invisible(x)
}


## Bruno Nicenboim (Tilburg, the Netherlands),
## Daniel Schad (Potsdam, Germany),
## Shravan Vasishth (Potsdam, Germany)
