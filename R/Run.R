

# library(shiny)
# library(bslib)
# library(DT)


# Upload large files:
# options(shiny.maxRequestSize = 30 * 1024^2)

#' @export
runApp = function(x = NULL) {
	# path = "...";
	# shiny::runApp(appDir = path);
	shiny::shinyApp(ui = getUI(), server = getServer())
}

