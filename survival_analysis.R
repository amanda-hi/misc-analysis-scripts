# Amanda Leonti
# October 12, 2018
# Purpose: Perform survival analysis using survival/survminer packages and create Kaplan-Meier plots.

#################################################################################

plot_survCurves <- function(data, test_type, event_col, censor_id = NULL, time_type, time_col, grouping_col = 1, group = NULL, palette = NULL, n_estimate = 5){
  
  #----------- Function arguments --------------------------------------------------------------------------#
  # data = Data frame containing survival information
  # test_type = A character string indicating the type of survival data being used (acceptable entries are OS, EFS, DFS, or RR, case insensitive)
  # event_col = A character string of the column name in the data that contains the event classifications
  # censor_id = A character (or numeric) string indicating the label used for censored data (i.e "Censored", 0, "censor", etc.)
  # time_type = A character string indicating the time format (options are "days" or "years")
  # grouping_col = A character string of the column name (in 'data') that contains the experimental groups.
  #                The default value of 1 is based on the survival package docs (chapter 2.1).
  # group = Plot title and/or specific subgroup that the curves are restricted to. Default is NULL
  # palette = Color palette to use for the plot, default is the "Paired" palette from Colorbrewer. Other options include:
  #                                                     - a vector of color names/RColorbrewer palette,
  #                                                     - "npg", "aaas", "lancet", "jco", "uchicago",
  #                                                       "ucscgb", "simpsons", "rickandmorty"
  # n_estimate = Single number (or vector of numbers) indicating the desired n-year survival estimate (ex. 5 for 5-year EFS estimate)
  #---------------------------------------------------------------------------------------------------------#
  
  
  
  #----------- Event coding info ---------------------------------------------------------------------------#
  # For OS events: The status indicator is normally 0 = alive, 1 = dead (but this only applies to OS events)
  # For EFS events/interval censored data: The status indicator is 0 = right censored, (fyi, the TARGET data is right censored)
  #                                                                1 = event at time, 
  #                                                                2 = left censored, 
  #                                                                3 = interval censored
  # For DFS events: Time to first event from EOI1 for patient in CR.
  #                 Event can be either "relapse" or "death w/o relapse", otherwise patient = censored.
  # For RR events: Similar to DFS, but calculated using the Competing Risks method (1-KM estimate). see here for more info:
  #                https://www.publichealth.columbia.edu/research/population-health-methods/competing-risk-analysis
  #                First events = either "primary" events (ex. relapse), or "secondary" events (ex. death w/o relapse)
  #                Non events = censored
  #                NOTE - level names cannot contain spaces or it will mess up the plot & curves
  #---------------------------------------------------------------------------------------------------------#
  
  library(cmprsk)
  library(dplyr)
  library(GGally)
  library(survival)
  library(survminer)
  
  ###########################################################################
  #----------------------- Functions ---------------------------------------#
  ###########################################################################
  
  #------------- Kaplan-Meier estimate plot ------------------#
  
  generate_KMplot <- function (survObj, test_type) {
    
    # Fitting the survival data to create Kaplan-Meier curves 
    # Implemented this fix for using global variables: https://github.com/kassambara/survminer/issues/342
    formula <- as.formula(paste0("survObj ~ ", grouping_col))
    fit <- survminer::surv_fit(formula = formula, data = data) # Using the formula generated above to fit the Kaplan-Meier curves
    surv_medians <- survminer::surv_median(fit) # Calculating median survival estimates & cleaning up strata labels
    surv_medians$strata <- gsub(".*=", "", surv_medians$strata)
    surv_medians$median <- ifelse(is.na(surv_medians$median), "N.R.", round(surv_medians$median, 2))
    surv_medians$strata_label <- ifelse(surv_medians$median == "N.R.", "N.R.", paste(surv_medians$median, time_type))
    
    # Creating the plot title using the group & test type information
    if (is.null(group)) {
      title = NULL # This isn't necessary, I've updated the x-axis labels to include this info
    } else {
      title = group
      # title = paste0("Kaplan-Meier curves for ", test_type, " (", group, ")")
    }
    
    # Fixing the names of the strata to make them more concise & descriptive
    if (grouping_col != 1) {
      names(fit$strata) <- gsub(".+=", "", names(fit$strata))
    } else {
      fit$strata <- fit$n
      names(fit$strata) <- "All"
    }
    
    # Adding patient count (n) to the strata name
    for (x in seq(length(names(fit$strata)))) {
      names(fit$strata)[x] <- paste0(names(fit$strata[x]), 
                                     " (n = ", fit$n[[x]], 
                                     ", median = ", surv_medians$strata_label[x], ")")
    }
    
    ylab <- case_when(test_type == "OS" ~ "Overall survival",
                      test_type == "EFS" ~ "Event-free survival",
                      test_type == "DFS" ~ "Disease-free survival",
                      test_type == "RR" ~ "Relapse risk")
    
    # Creating the Kaplan-Meier plot
    plot <- ggsurvplot(fit, data = data,
                       palette = palette,
                       pval = ifelse(grouping_col == 1, FALSE, TRUE),
                       pval.coord = c(0.2, 0.08),
                       ggtheme = theme_classic() + 
                         theme(plot.title = element_text(hjust = 0.5, size = 12),
                               axis.text = element_text(size = 11),
                               axis.title = element_text(size = 12), 
                               legend.text = element_text(size = 8)),
                       legend = "bottom",
                       title = title, 
                       surv.scale = "percent") +
      guides(fill = guide_legend(title = NULL, ncol = legendDim),
             color = guide_legend(title = NULL, ncol = legendDim)) +
      labs(x = paste0("Time (", time_type, ")"), y = ylab)
    
    plot
    
    return(list(plot = plot$plot, fit = fit))
  }
  
  
  #------------- Cumulative incidince of risk plot -----------#
  
  generate_CRplot <- function(survObj) {
    
    # Creating a Competing Risk plot
    plot <- ggcompetingrisks(fit = survObj, 
                             multiple_panels = F,
                             palette = palette,
                             legend = "bottom", 
                             ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 12))) +
      guides(fill = guide_legend(title = "Event"),
             color = guide_legend(title = "Event")) +
      labs(x = paste0("Time (", time_type, ")"), title = "Cumulative Incidence Estimate for Relapse Risk (RR)") +
      scale_linetype_manual(name = "Group", values = rep(c("solid", "dashed"), times = nstrata/2))
    
    plot
  }
  
  #----------- Cox PH Plot --------------------------#
  
  coxPH_fit <- function(survObj, grouping_col){
    
    # Running a Cox Proportional Hazards regression
    formula <- as.formula(paste0("survObj ~ ", grouping_col))
    coxPH_model <- coxph(formula, data = data)
    
    # Visualizing the results using a forest plot
    # source("Modified_ggforest_Plot.R")
    # forestPlot <- ggforest2(coxPH_model)
    
    return(list(Model = coxPH_model))
  }
  
  palette_picker <- function(x) {
    # See https://stackoverflow.com/questions/62005334/why-cant-case-when-return-different-length-vectors for reason why the list() is needed.
    p <- case_when(x >= 8 ~ list(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(x)),
                   x < 8 ~ list(RColorBrewer::brewer.pal(9, "Set1")[c(1:5, 7:9)]))
    p <- unlist(p)
    return(p)
  }
  
  
  ###########################################################################
  #----------------------- Script body -----------------------------------#
  ###########################################################################
  
  # Performing some data prep before plotting
  data <- as.data.frame(lapply(data, unlist), optional = T) # Ensuring that no dataframe columns/values are actually lists
  
  # Reconfigures the legend dimensions if there are a lot of groups
  nstrata <- ifelse(grouping_col != 1, length(unique(data[[grouping_col]])), 1)
  legendDim <- ifelse(nstrata > 10, floor(nstrata/2), 1)
  
  if (grouping_col != 1) {
    test_data <- as.factor(data[,grouping_col]) # Making sure this variable is a factor, otherwise Surv() won't work
    
    # Using number of individual group levels to select a color palette
    nlev <- length(levels(test_data))
    palette <- if (is.null(palette)) {
      palette_picker(nlev)
    } else {
      palette
    }
    
  } else {
    test_data <- NA
    palette <- if (is.null(palette)) {
      colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(1)
    } else {
      palette
    }
  }
  
  # Converting the time unit to years, instead of days
  if (time_type == "days") {
    time_col <- as.numeric(data[,time_col])/365
    time_type <- "years" # Updating variable
  } else {
    time_col <- data[,time_col]
  }
  
  test_type <- toupper(test_type) # Casting to uppercase so that user entry can be case insensitive
  
  if (test_type == "RR") {
    survObj <- cuminc(ftime = as.numeric(time_col), # Creating the Competing Risks survival object 
                      fstatus = data[,event_col],
                      cencode = censor_id, 
                      group = test_data)
    
    plot <- generate_CRplot(survObj)
    
  } else {
    survObj <- Surv(time = as.numeric(time_col),    # Creating the Kaplan-Meier survival object
                    event = as.numeric(data[,event_col]))
    
    # This test determines whether or not the groups differ significantly, in regard to survival
    if (grouping_col != 1) {
      logRank.formula <- as.formula(paste0("survObj ~ ", grouping_col))
      logRank.test <- survdiff(logRank.formula, data = data) 
    } else {
      logRank.test <- NA
    }
    
    # Creating the Kaplan-Meier plot using the Surv object from above
    plot <- generate_KMplot(survObj = survObj, 
                            test_type = test_type)
    
    # Fitting a Cox PH regression model (also using the Surv object from above)
    coxRegression <- coxPH_fit(survObj = survObj, grouping_col = grouping_col)
    
    n_estimate <- summary(plot$fit, time = n_estimate, extend = T)$surv
    names(n_estimate) <- names(plot$fit$strata)
  }
  
  
  # ---------- Returning final plots and model objects ------------------------------------------#
  
  if (test_type == "RR") {
    return(list(cr_plot = plot, cuminc_object = survObj))
  } else {
    return(list(km_plot = plot$plot, cox_model = coxRegression, nYear_estimate = n_estimate, fitted_curve = plot$fit, surv_object = survObj, log_rank = logRank.test))
  }
  
}
