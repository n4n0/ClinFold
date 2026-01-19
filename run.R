# ============================================================
# ClinFold Installation and Run Script
# ============================================================

cat("\n========================================\n")
cat("ClinFold: Structural Impact Viewer\n")
cat("for Genetic Variants\n")
cat("========================================\n\n")

required_packages <- c(
  "shiny", "httr", "jsonlite", "dplyr", "tibble",
  "DT", "bio3d", "purrr", "glue", "stringr", 
  "ggplot2", "plotly", "NGLVieweR"
)

cat("Checking required packages...\n\n")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat("  OK", pkg, "\n")
  }
}

cat("\n========================================\n")
cat("Starting ClinFold...\n")
cat("========================================\n\n")

shiny::runApp("app.R", launch.browser = TRUE)