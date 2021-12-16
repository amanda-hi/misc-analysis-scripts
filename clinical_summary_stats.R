# Amanda Leonti
# Purpose: Generate clinical characteristic summary tables similar to those provided by 
# COG statisticians for use in manuscripts, analyses, and presentations.
#############################################

# Additional summary stats packages to test:
# https://thatdatatho.com/easily-create-descriptive-summary-statistic-tables-r-studio/

generate_sumStats_table <- function(df, cols, group_var, title = NULL, write_to_file = FALSE, filename = NULL, rowWise = FALSE) {
  
  # df = Dataframe containing clinical data to be summarized
  # cols = Vector containing the names of columns/variables (as strings) to perform tests on & include in the final table
  # group_var = column name (as a string) of the categorical variable used to test
  
  library(arsenal) 
  library(openxlsx)
  library(dplyr)

  group_var <- sym(group_var)
  
  df <- df %>%
    dplyr::select(any_of(cols), any_of(group_var))
  
  if (rowWise == TRUE) {
    tests <- df %>%
      tableby(eval(group_var) ~ ., data = ., cat.stats = "countrowpct", # Performs rowwise calculations across selected characteristics
              total.pos = "before",                                     # and patient categories (specified by the "cols" variable), rather than within
              test.pname = "p-value",                                   # patient groups (specified by the "group_var" variable),
              cat.test = "chisq",                                       # see https://github.com/mayoverse/arsenal/issues/143#issuecomment-427933192
              digits = 2,
              cat.simplify = FALSE, 
              ordered.simplify = FALSE,
              numeric.stats = c("Nmiss", "meansd", "median", "range"))
  } else {
    tests <- df %>%
      tableby(eval(group_var) ~ ., data = .,                          # Performs analysis in a column-wise fashion (aka the default) 
              total.pos = "before",                                   # within specified patient groups (the "group_var") variable
              test.pname = "p-value", 
              cat.test = "chisq", 
              digits = 2,
              cat.simplify = FALSE, 
              ordered.simplify = FALSE,
              numeric.stats = c("Nmiss", "meansd", "median", "range"))
  }
  
  table <- as.data.frame(summary(tests, title = "Clinical Summary", text = TRUE))
  
  
  
  ####### OPTIONAL: Writing table out to an Excel file #########################
  
  if (write_to_file == TRUE) {
    
    if (is.null(filename)) {
      stop("Please specify a filename! Include the full path to the destination directory, otherwise it will be written to the current working directory.")
    }
    
    write.xlsx(table, file = filename, 
               firstRow = TRUE, 
               borders = "columns",
               # Column header styling
               headerStyle = createStyle(halign = "center", textDecoration = "Bold", border = "bottom"),
               overwrite = T)
  }
  
  return(table)
}