# Amanda Leonti
# Nov 13th, 2019
# Purpose: the "ggforest" function from the survminer package is not very customizable, 
# so I grabbed the code and modified it a bit to make it my own. 

########################################################################

library(broom)
library(grid)

# Function defined by the package developers to pull data out of a fitted Cox PH model object
.get_data <- function(fit, data = NULL, complain = TRUE) {
  if(is.null(data)){
    if (complain)
      warning ("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
    if (is.null(data))
      stop("The `data` argument should be provided either to ggsurvfit or survfit.")
  }
  data
}

ggforest2 <- function(model, data = NULL,
                      main = "Hazard ratio", 
                      cpositions = c(0.02, 0.22, 0.4),
                      continuous = FALSE,
                      fontsize = 0.7, refLabel = "reference", 
                      noDigits = 2, boxColor = "black", sigColor = "black") {

  
  ####################### DATA SETUP #############################
  
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(class(model) == "coxph")
  
  # get data and variables/terms from cox model
  data  <- .get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  
  # use broom to get some required statistics
  coef <- as.data.frame(broom::tidy(model))
  gmodel <- glance(model)
  
  # extract statistics for every variable
  allTerms <- lapply(seq_along(terms), function(i){
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data),
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data),
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[,1:2], 1, paste0, collapse="")
  
  # use broom again to get remaining required statistics
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds,])[,c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high", "pos")]
  toShowExp <- toShow[,5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits=noDigits)
  toShowExpClean <- data.frame(toShow,
                               pvalue = signif(toShow[,4],noDigits+1),
                               toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, noDigits+1), " ",
                                 ifelse(toShowExpClean$p.value < 0.05, "*",""),
                                 ifelse(toShowExpClean$p.value < 0.01, "*",""),
                                 ifelse(toShowExpClean$p.value < 0.001, "*",""))
  toShowExpClean$ci <- paste0("(",toShowExpClean[,"conf.low.1"]," - ",toShowExpClean[,"conf.high.1"],")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  
  # make label strings
  toShowExpClean$N <- paste0("(N=",toShowExpClean$N,")")
  
  # flip order
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  
  # make plot twice as wide as needed to create space for annotations
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  
  # increase white space on right for p-vals
  rangeplot[2] <- rangeplot[2] + .35 * diff(rangeb)
  width <- diff(rangeplot)
  
  # y-coordinates for label
  y_variable <- rangeplot[1] +  cpositions[1] * width
  
  if (continuous == TRUE) {
    y_nlevel <- rangeplot[1]  +  cpositions[2] * width
  } else {
    y_nlevel <- 0
    y_cistring <- rangeplot[1]  +  cpositions[2] * width
  }
  
  y_stars <- rangeplot[2]/1.5
  x_annotate <- seq_len(nrow(toShowExpClean))
  
  # the units for geom_text fontsize are mm (https://github.com/tidyverse/ggplot2/issues/1828)
  annot_size_mm <- fontsize *
    as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))
  
  
  
  ####################### PLOT GENERATION #############################
  
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
    
    ggtitle(main) +
    xlab("") +
    theme_light() +
    
    geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                  ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                  fill = ordered(seq_along(var) %% 2 + 1))) +
                # fill = var)) +
    
    # Dotted line to mark an HR of 1
    geom_hline(yintercept = 1, linetype = 3, colour = "grey35") +
    
    # Colors for the alternating grey/white background bars
    scale_fill_manual(values = c("#FFFFFF33", "grey95"), guide = "none") +
    
    # Boxes for the forest plot that represent HR values
    geom_point(pch = 15, size = 4, color = boxColor) +
    
    # CI intervals
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
    
    # Flips the plot sideways
    coord_flip(ylim = exp(rangeplot)) +
    
    scale_y_log10(
      name = "",
      labels = sprintf("%g", breaks),
      expand = c(0.02, 0.02),
      breaks = breaks) +

    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = rel(1)),
          plot.title = element_text(hjust = 0.5, size = annot_size_mm*3, face = "bold")) +

    # Vertical columns of HR, CI, and p-value data
    annotate(geom = "text", x = x_annotate, y = exp(y_variable),
             label = toShowExpClean$var, fontface = "bold", hjust = 0,
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0,
             label = toShowExpClean$level, vjust = -0.1, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShowExpClean$estimate.1, size = annot_size_mm,
             vjust = 0.5) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShowExpClean$ci, hjust = -0.4, vjust = 0.5, size = annot_size_mm/1.4,  fontface = "italic") +
    # annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
    #          label = toShowExpClean$ci, size = annot_size_mm/1.4,
    #          vjust = 1.1,  fontface = "italic") +
    annotate(geom = "text", x = x_annotate, y = exp(y_stars),
             label = toShowExpClean$stars, size = annot_size_mm,
             hjust = -0.2,  fontface = ifelse(toShowExpClean$stars < 0.05, "bold.italic", "plain"), 
             colour = ifelse(toShowExpClean$stars < 0.05, sigColor, "black")) +
    annotate(geom = "text", x = 0.5, y = exp(y_variable),
             label = paste0("Global p-value (Log-Rank): ", format.pval(gmodel$p.value.log, eps = ".001"), ", AIC: ", round(gmodel$AIC,2)),
             size = annot_size_mm/1.4, hjust = 0, vjust = 1.2,  fontface = "italic")
  # annotate(geom = "text", x = 0.5, y = exp(y_variable),
  #          label = paste0("# Events: ", gmodel$nevent, "; Global p-value (Log-Rank): ",
  #                         format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", round(gmodel$AIC,2),
  #                         "; Concordance Index: ", round(gmodel$concordance,2)),
  #          size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")
}