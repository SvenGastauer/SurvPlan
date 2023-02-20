#' Run the Survey Palnner app
#' @import shiny
#' @import leaflet
#' @import DT
#' @import ggplot2
#' @import sf
#' @import purrr
#' @import dplyr
#' @import sp
#' @import tidyr
#' @export
#'
planner<- function(){
  appDir <- system.file("planner", package = "SurvPlan")
  if (appDir == "") {
    stop("Could not find directory. Try re-installing `SurvPlan`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal",launch.browser = TRUE)
}

