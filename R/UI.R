

### UI Functions

init.theme = function() {
    bs_theme(
        bg = "#101010", 
        fg = "#FDF7F7", 
        primary = "#ED79F9"
    );
}


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
	### Shiny functions
	# - NO need to be visible in the app;
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
	# <----------------->
	### UI:
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
						value = 0.2, min = 0, max = 5, step = 0.25),
					selectInput("fltAllele",
						label = "Filter alleles:",
						choices = list("All" = "All", "Common / Known" = "Common",
							">= 3%" = ">= 3%", ">= 1%" = ">= 1%",
							"Uncommon (< 1%)" = "< 1%", "Uncommon (NA)" = "Uncommon"),
						selected = "All"),
					checkboxInput("chkRegex", "Regex Search: Data", value = TRUE),
					downloadButton("downloadData", "Download"),
				),
				mainPanel(DT::DTOutput("tblData"))
		)),
		tabPanel("Alleles", # icon = icon("HLA"),
			mainPanel(DT::DTOutput("tblAlleles"))
		),
		# Population Coverage
		tabPanel("Epitopes", # icon = icon("Epitopes"),
			fluidRow(
			column(4,
				actionButton("btnCovHLA", "Display HLA"),
				textOutput("txtBtnDisplay")),
			column(8, DT::DTOutput("tblTotalPopulation"))
			),
			fluidRow(br()),
			fluidRow(
			column(9, DT::DTOutput("tblPeptides")),
			column(3, DT::DTOutput("tblAllelesPP"))
			),
			fluidRow(
				column(5,
				downloadButton("downloadPP", "Download"),
				actionButton("printPPSel", "Print Selection"),
				actionButton("btnRemainingEpi", "Remaining"),
				),
				column(7,
				fluidRow("Ti = Population coverage (assuming independence of A/B/C-alleles)"),
				fluidRow("Tn = Simple Total (naive sum)")
				)
			),
			fluidRow(
			column(12, DT::DTOutput("tblRemainingEpi"))
			)
		),
		# Sub-Sequences
		tabPanel("SubSeq", # icon = icon("SubSeq"),
			fluidRow(
			column(4,
				textInput("inSubSeq", "Peptide Seq", value = "")
			),
			column(4,
				actionButton("btnSearchSubSeq", "Search"),
				textOutput("txtBtnSearchSubSeq"))
			),
			fluidRow(
			column(12, DT::DTOutput("tblSummarySubSeq"))
			),
			fluidRow(
			column(12, DT::DTOutput("tblSubSeq"))
			)
		),
		tabPanel("EpiSummary", # icon = icon("Summary"),
			fluidRow(
			column(4, # Epitopes which will be summarized
				textInput("inEpiSummary", "Epitopes", value = "")
			),
			column(4,
				actionButton("btnEpiSummary", "Summary"),
				textOutput("txtBtnEpiSummary"))
			),
			fluidRow(
			column(12, DT::DTOutput("tblEpiSummary"))
			)
		),
		tabPanel("Protein", # icon = icon("Protein"),
				plotOutput("imgProtein")
		),
		tabPanel("Help", # icon = icon("Help"),
				uiOutput("txtHelp")
		)
	)))
}
