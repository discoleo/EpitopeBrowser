

### UI Functions

init.theme = function() {
    bs_theme(
        bg = "#101010", 
        fg = "#FDF7F7", 
        primary = "#ED79F9"
    );
}

### Shiny functions
# Layout:
sidebarLayout = function(...) shiny::sidebarLayout(...);
fluidRow = function(...) shiny::fluidRow(...);
column = function(...) shiny::column(...);
# Panels:
tabPanel = function(...) shiny::tabPanel(...);
mainPanel = function(...) shiny::mainPanel(...);
titlePanel = function(...) shiny::titlePanel(...);
sidebarPanel = function(...) shiny::sidebarPanel(...);
# Input:
sliderInput = function(...) shiny::sliderInput(...);
selectInput = function(...) shiny::selectInput(...);
checkboxInput = function(...) shiny::checkboxInput(...);

fileInput.csv = function(id, label) {
	shiny::fileInput(id, label,
		multiple = FALSE,
		accept = c(
			"text/csv",
			"text/comma-separated-values,text/plain",
			".csv", ".txt")
    );
}

getUI = function() {
	shiny::shinyUI(
	shiny::fluidPage(
	shiny::navbarPage("HLA-Browser", id="menu.top",
		# theme = init.theme(),
		tabPanel("Data", # icon = icon("upload file"),
			sidebarLayout(
				sidebarPanel(
					titlePanel("Load file"),
					fileInput.csv("file", "Select CSV file"),
					sliderInput(inputId = "rank", label = "Rank",
						value = 0.2, min = 0, max = 2, step = 0.025),
					selectInput("fltAllele",
						label = "Filter alleles:",
						choices = list("All" = "all", "Common" = "common",
							"Uncommon" = "uncommon"),
						selected = "all"),
					checkboxInput("chkRegex", "Regex Search: Data", value = TRUE)
				),
				mainPanel(DT::DTOutput("tblData"))
		)),
		tabPanel("Alleles", # icon = icon("HLA"),
			mainPanel(DT::DTOutput("tblAlleles"))
		),
		tabPanel("Epitopes", # icon = icon("Epitopes"),
			fluidRow(
			column(4, actionButton("ppHLA", "Display HLA")),
			column(8, DT::DTOutput("tblTotalPopulation"))
			),
			fluidRow(br()),
			fluidRow(
			column(9, DT::DTOutput("tblPeptides")),
			column(3, DT::DTOutput("tblAllelesPP"))
			),
			fluidRow(
				downloadButton("downloadPP", "Download")
			)
		)
	)
	))
}
