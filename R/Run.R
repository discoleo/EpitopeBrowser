

# library(shiny)
# library(bslib)
# library(DT)

#' @export
runApp = function(x = NULL) {
	# path = "...";
	# shiny::runApp(appDir = path);
	shiny::shinyApp(ui = getUI(), server = getServer())
}

