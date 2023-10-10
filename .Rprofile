source("renv/activate.R")
if (interactive() && requireNamespace("rsthemes", quietly = TRUE)) {
  # Set preferred themes if not handled elsewhere..
  rsthemes::set_theme_light("One Light {rsthemes}")  # light theme
  rsthemes::set_theme_dark("One Dark {rsthemes}") # dark theme
  
  # Whenever the R session restarts inside RStudio...
  setHook("rstudio.sessionInit", function(isNewSession) {
    # Automatically choose the correct theme based on time of day
    rsthemes::use_theme_auto()
  }, action = "append")
}
