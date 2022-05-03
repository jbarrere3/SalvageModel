# Install targets package if needed
if(!("targets" %in% installed.packages())) install.packages("targets")
# Load targets
library(targets)
# Launch the targets pipeline
tar_make()
