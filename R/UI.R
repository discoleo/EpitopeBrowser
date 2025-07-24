

### UI Functions

init.theme = function() {
    bs_theme(
        bg = "#101010", 
        fg = "#FDF7F7", 
        primary = "#ED79F9"
    );
}



# version = 1: Old variant;
getUI = function(versionPop = 2) {
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
					sliderInput(inputId = "fltRank", label = "Rank",
						value = 0.2, min = 0, max = 5, step = 0.25),
					checkboxInput("chkRegex", "Regex Search: Data", value = TRUE),
					downloadButton("downloadData", "Download"),
					actionButton("btnDataStats", "Stats"),
					# Filter: AA Position
					fluidRow(
					column(6, textInput("fltSeqStart", "Start", "", width = 150)),
					column(6, textInput("fltSeqEnd", "End", "", width = 150))
					),
					# HLA Alleles: Region
					inputHLARegion(),
					inputHLAFreq(),
					fluidRow("HLA: Frequency-data extracted from",
						"various Regions of interest."
					),
				),
				mainPanel(
				fluidRow(DT::DTOutput("tblData")),
				# Stats:
				fluidRow(htmlOutput("txtDataStats")),
				fluidRow(DT::DTOutput("tblDataStats")),
				)
		)),
		### HLA Alleles
		tabPanel("Alleles", # icon = icon("HLA"),
			mainPanel(DT::DTOutput("tblAlleles"))
		),
		### Population Coverage
		# - see file: UI.PopCover.R;
		getUI.PopCover(version = versionPop),
		### Sub-Sequences
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
			fluidRow(DT::DTOutput("tblEpiSummaryUnique")),
			fluidRow(
			column(12, DT::DTOutput("tblEpiSummary"))
			),
		),
		# Regions: Compare HLA Freq
		tabPanel("Regions", # icon = icon("Regions"),
			fluidRow(DT::DTOutput("tblRegionsHLA"))
		),
		# View Epitopes on Protein
		tabPanel("Protein", # icon = icon("Protein"),
			plotOutput("imgProtein")
		),
		tabPanel("Help", # icon = icon("Help"),
			uiOutput("txtHelp")
		)
	)))
}

#######################

# Input: File Open
fileInput.csv = function(id, label) {
	shiny::fileInput(id, label,
		multiple = FALSE,
		accept = c(
			"text/csv",
			"text/comma-separated-values,text/plain",
			".csv", ".txt")
    );
}

### Input: Allele Region
inputHLARegion = function() {
	selectInput("fltAlleleRegion",
		label = "HLA Set:",
		choices = list(
			"Germany" = "De", "Italy" = "Italy",
			"Mix" = "Mix", "Hungary" = "Hungary",
			"Ro (Fundeni)" = "Ro-Fundeni",
			"Not yet" = "Not yet"),
		selected = "Germany");
}

### Input: Allele Freq-Filter
inputHLAFreq = function() {
	selectInput("fltAllele",
		label = "Filter alleles:",
		choices = list("All" = "All", "Common / Known" = "Common",
					">= 3%" = ">= 3%", ">= 1%" = ">= 1%",
					"Uncommon (< 1%)" = "< 1%", "Uncommon (NA)" = "Uncommon"),
		selected = "All");
}
