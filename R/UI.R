

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
				mainPanel(DT::DTOutput("tblData"))
		)),
		tabPanel("Alleles", # icon = icon("HLA"),
			mainPanel(DT::DTOutput("tblAlleles"))
		),
		# Population Coverage
		getUI.PopCover(version = versionPop),
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

#######################

### Population Coverage

getUI.PopCover = function(version = 2) {
	if(version == 2) {
		return(getUI.PopCover.v2());
	} else if(version == 1) {
		return(getUI.PopCover.v1());
	} else if(version == 0) {
		return(getUI.PopCover.old());
	}
	warning("Unsupported version!");
}
getUI.PopCover.v2 = function() {
	tabPanel("Epitopes", # icon = icon("Epitopes"),
	fluidRow(
		column(4,
			fluidRow("Ti = Population coverage (assuming independence of A/B/C-alleles)"),
			fluidRow("Tn = Simple Total (naive sum)"),
			fluidRow("Ta = All alleles included"),
		),
		column(8, DT::DTOutput("tblTotalPopulation"))
	),
	fluidRow(br()),
	fluidRow(
		column(9,
			fluidRow(DT::DTOutput("tblPeptides"))
		),
		# Right "Panel": less mouse movement
		column(3,
		fluidRow(
			column(1, ""),
			column(11, offset = 1,
			fluidRow(
				downloadButton("btnDownloadPP", "Download"),
				actionButton("printPPSel", "Print Selection"),
			),
			fluidRow(h3("Population Coverage:")),
			fluidRow(
				textOutput("txtBtnDisplay"),
				actionButton("btnCovHLA", "Display HLA"),
			),
			fluidRow(h3("Remaining Epitopes:")),
			fluidRow(
				actionButton("btnRemainingEpi", "Remaining"),
				actionButton("btnGoToPage", "Go To"),
			),
			fluidRow(textOutput("txtPPTblPage"))
			))
		)),
	fluidRow(
		column(3, DT::DTOutput("tblAllelesPP")),
		column(9, DT::DTOutput("tblRemainingEpi"))
	)
	)
}

### Variant:
# - a lot of mouse movement;
getUI.PopCover.v1 = function() {
	tabPanel("Epitopes", # icon = icon("Epitopes"),
	fluidRow(
		column(4,
			fluidRow("Ti = Population coverage (assuming independence of A/B/C-alleles)"),
			fluidRow("Tn = Simple Total (naive sum)"),
			fluidRow("Ta = All alleles included"),
		),
		column(8, DT::DTOutput("tblTotalPopulation"))
	),
	fluidRow(br()),
	fluidRow(
		# Left "Panel"
		column(3,
			fluidRow(
				downloadButton("btnDownloadPP", "Download"),
				actionButton("printPPSel", "Print Selection"),
			),
			fluidRow(h3("Population Coverage:")),
			fluidRow(
				textOutput("txtBtnDisplay"),
				actionButton("btnCovHLA", "Display HLA"),
			),
			fluidRow(h3("Remaining Epitopes:")),
			fluidRow(
				actionButton("btnRemainingEpi", "Remaining"),
				actionButton("btnGoToPage", "Go To"),
			),
			fluidRow(textOutput("txtPPTblPage")),
		),
		column(9,
			fluidRow(DT::DTOutput("tblPeptides"))
		)),
	fluidRow(
		column(9, DT::DTOutput("tblRemainingEpi")),
		column(3, DT::DTOutput("tblAllelesPP"))
	)
	)
}

# Old variant:
getUI.PopCover.old = function() {
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
			fluidRow(
				downloadButton("btnDownloadPP", "Download"),
				actionButton("printPPSel", "Print Selection"),
				actionButton("btnRemainingEpi", "Remaining"),
				actionButton("btnGoToPage", "Go To"),
			),
			fluidRow(textOutput("txtPPTblPage"))
		),
		column(7,
			fluidRow("Ti = Population coverage (assuming independence of A/B/C-alleles)"),
			fluidRow("Tn = Simple Total (naive sum)")
		)
	),
	fluidRow(
		column(12, DT::DTOutput("tblRemainingEpi"))
	)
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
