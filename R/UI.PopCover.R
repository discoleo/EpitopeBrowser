
### Population Coverage
# UI: Tab for Pop. Cover;

# Select Variant:
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

### Variant v2
# - improved UX: Panel on right;
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
				actionButton("btnCovClear", "Clear"),
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
				actionButton("btnCovClear", "Clear"),
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

