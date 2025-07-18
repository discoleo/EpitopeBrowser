

getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		fltRank   = 0.55, # Default value for Rank-Filter;
		hla.strip = TRUE, # HLA-A... => A...;
		sep = ",",        # csv Separator
		HLA = hla(type = "De"),
		HLA.trim.2A = TRUE, # Trim: D[PQR]A*...;
		# Regex & Other Options:
		reg.Data  = TRUE,  # Regex for Data-Table
		reg.PP    = TRUE,  # Regex for Epitopes-Table
		highlight = TRUE,  # Highlight Search Term
		# Sub-Seq:
		allEpi.SubSeq = TRUE,
		allEpi.EpiSummary = TRUE,
		# Protein Graph:
		col.Pr    = "#FF0032A0",
		border.Pr = "#640000A0",
		lwd.Pr   = 1.5 # NOT used with Polygon
	);
	
	const = list(
		DisplayHLA     = list(
			Warn     = "Select first some Epitopes in the table below.",
			More     = "Additional epitopes can be selected."),
		SearchSubSeq   = "Search epitopes within this peptide.",
		SearchNoSubSeq = paste("No epitopes found.",
			"Please enter another peptide sequence for a new search."),
		EpiSummary     = "Enter list of epitopes to summarize.",
		NoEpiSummary   = paste("No epitopes found.",
			"Please enter a new list of epitopes to summarize."),
		NULL
	);
	
	### Init:
	updateNumericInput(session, "rank", value = options$fltRank);
	
	# Dynamic variable
	values = reactiveValues(
		Active    = "Not",  # Active Tab
		fullData  = NULL,   # initial Data
		dfGlData  = NULL,   # Globally filtered Data
		dfFltData = NULL,   # Data filtered in Table
		reg.Data  = options$reg.Data,
		fltRank   = NULL,   # is set automatically
		fltAllele = NULL,
		fltCols   = NULL,   # TODO
		fltSeqPos = NULL,   # AA-Position in Seq
		multSeq   = FALSE,
		typeHLA   = 1,      # HLA Class
		# Population Coverage:
		dfPopCoverPP = NULL,
		dfPopAlleles = NULL, # Alleles Covered by selected Set
		dfTotalPopulation = NULL, # HLA & Epitope
		fltHLAEpiSel      = NULL, # Only HLA covered by Selection
		dfRemainingEpi    = NULL,
		optRemainingEpi   = list(), # Options/Filters for the table
		pageTblPP    = NULL,      # Page in main PopCoverage Table
		# Sub-Sequences:
		warnSubSeq = FALSE,
		warnEpiSummary = FALSE,
		NULLARG = NULL
	);
	
	
	### Menu Tabs
	observeEvent(input$menu.top, {
		if(input$menu.top != "Data") {
			print("Switched");
			filter.byTable();
			# values$Active = "Other";
		}
	})
	
	# Reset Filters on Data table
	hasData = function() { ! is.null(values$dfGlData); }
	
	### Population Coverage: each Epitope
	freq.hla = function() {
		x = values$dfFltData[c("Peptide", "HLA")];
		freq.population(x, options$HLA, type = values$typeHLA);
	}
	
	trim = function(x) {
		sub("(?i)^HLA-", "", x);
	}
	filter.df = function() {
		if(is.null(values$fullData)) return();
		getFilter.tblData();
		lim.rank  = values$fltRank;
		fltAllele = values$fltAllele;
		x = values$fullData;
		# Trick: lim = 0 => accept ALL Ranks;
		if(lim.rank > 0) x = x[x$Rank <= lim.rank, ];
		x = filter.HLA(x, fltAllele);
		if(! is.null(values$fltSeqPos)) x = filter.seqPos(x);
		#
		values$dfGlData = x;
		cat("Rows: ", nrow(x), "\n");
	}
	filter.seqPos = function(x) {
		flt  = values$fltSeqPos;
		fltS = flt$Start;
		fltE = flt$End;
		isS = ! is.null(fltS);
		isE = ! is.null(fltE);
		isSeq = rep(TRUE, nrow(x));
		# Note: Initial names;
		if(isS) {
			isSeq = (x$start >= fltS[1]);
			if(length(fltS) > 1)
				isSeq = isSeq & (x$start <= fltS[2]);
		}
		if(isE) {
			isLim1 = length(fltE) == 1;
			isSeq  = isSeq &
				if(isLim1) (x$end <= fltE[1])
				else ((x$end >= fltE[1]) & (x$end <= fltE[2]));
		}
		x = x[isSeq, ];
		return(x);
	}
	filter.HLA = function(x, fltAllele) {
		if(fltAllele == "All") return(x);
		#
		hla = if(options$hla.strip) x$HLA else trim(x$HLA);
		if(fltAllele == ">= 3%") {
			frequentHLA = options$HLA$HLA[options$HLA$Freq >= 0.03];
			isHLA = hla %in% frequentHLA;
		} else if(fltAllele == ">= 1%") {
			frequentHLA = options$HLA$HLA[options$HLA$Freq >= 0.01];
			isHLA = hla %in% frequentHLA;
		} else if(fltAllele == "< 1%") {
			frequentHLA = options$HLA$HLA[options$HLA$Freq < 0.01];
			isHLA = hla %in% frequentHLA;
		} else {
			isHLA = hla %in% options$HLA$HLA;
			if(fltAllele == "Uncommon") isHLA = ! isHLA;
		}
		# print("HLA"); print(head(options$HLA))
		x = x[isHLA, ];
		return(x);
	}
	# Table Filter:
	filter.byTable0 = function() {
		id = input$tblData_rows_all;
		if(is.null(id)) return(NULL);
		values$dfGlData[id, ];
	}
	filter.byTable = function() {
		# Analysis is performed on the filtered data;
		print("Filter Table 2");
		id = input$tblData_rows_all;
		values$dfFltData = values$dfGlData[id, ];
	}
	#
	option.regex = function(x, varia = NULL, caseInsens = TRUE) {
		opt = list(search = list(regex = x, caseInsensitive = caseInsens),
			searchHighlight = options$highlight);
		if( ! is.null(varia)) opt = c(opt, varia);
		return(opt);
	}
	
	# File Input
	observeEvent(input$file, {
		file1 = input$file;
		if (is.null(file1))
			return(NULL);
		#
		x = read.epi(file1$datapath,
			hla.strip = options$hla.strip, sep = options$sep);
		# Multiple Protein Sequences:
		values$multSeq = if("Seq" %in% names(x)) TRUE else FALSE;
		# HLA Class 2:
		if(any(grepl("(?i)^D", x$HLA))) {
			values$typeHLA = 2;
		} else values$typeHLA = 1;
		# NO population data
		if(options$HLA.trim.2A) {
			x$HLA = trim.hla2A(x$HLA);
		}
		#
		values$fullData = x;
	})
	
	### Options: Data
	observeEvent(values$fullData, {
		filter.df();
	})
	observeEvent(input$rank, {
		if(! is.null(values$fltRank) && input$rank == values$fltRank) return();
		values$fltRank = input$rank;
		filter.df();
	})
	observeEvent(input$fltAllele, {
		if(! is.null(values$fltAllele) && input$fltAllele == values$fltAllele) return();
		values$fltAllele = input$fltAllele;
		filter.df();
	})
	observeEvent(input$chkRegex, {
		isReg = input$chkRegex;
		if(values$reg.Data != isReg) {
			values$reg.Data = isReg;
		}
	})
	
	### Flt: AA-Position
	observeEvent(input$fltSeqStart, {
		txt = input$fltSeqStart;
		if(nchar(txt) == 0) {
			if(is.null(values$fltSeqPos)) return();
			limE = values$fltSeqPos$End;
			if(is.null(limE)) {
				values$fltSeqPos = NULL;
			} else {
				values$fltSeqPos$Start = NULL;
			}
			filter.df();
			return();
		}
		# Values:
		val = splitRange(txt);
		if(length(val) == 0) {
			cat("Invalid Start position!\n");
			return();
		}
		values$fltSeqPos$Start = val;
		filter.df();
	})
	observeEvent(input$fltSeqEnd, {
		txt = input$fltSeqEnd;
		if(nchar(txt) == 0) {
			if(is.null(values$fltSeqPos)) return();
			limS = values$fltSeqPos$Start;
			if(is.null(limS)) {
				values$fltSeqPos = NULL;
			} else {
				values$fltSeqPos$End = NULL;
			}
			filter.df();
			return();
		}
		# Values:
		val = splitRange(txt);
		if(length(val) == 0) {
			cat("Invalid End position!\n");
			return();
		}
		values$fltSeqPos$End = val;
		filter.df();
	})
	
	### Tbl Filter
	getFilter.tblData = function() {
		flt = input$tblData_search_columns;
		flt = c("", flt); # Row ID
		if(all(flt == "")) {
			values$fltCols = NULL;
			return(NULL);
		}
		flt = lapply(flt, function(x) if(x == "") NULL else list(search = x));
		values$fltCols = flt;
		return(flt);
	}
	
	### Tables
	
	# Data
	dataTable = function() {
		# print("Rendering table!");
		flt = values$fltCols; # getFilter.tblData();
		if(! is.null(flt)) flt = list(searchCols = flt);
		if(is.null(values$dfGlData)) return();
		DT::datatable(values$dfGlData, filter = 'top',
			options = option.regex(values$reg.Data, varia = flt)) |>
		formatRound(c('Score', 'Rank'), c(4,2));
	}
	
	output$tblData = DT::renderDT(dataTable())
	
	# HLA Alleles
	output$tblAlleles <- DT::renderDT ({
		# values$Active = "Alleles";
		x = values$dfFltData$HLA;
		x = data.frame(table(x));
		names(x)[1] = "HLA";
		x = merge.hla(x, options$HLA);
		DT::datatable(x, filter = 'top') |>
			formatRound('Population (%)', 2);
	})
	
	### Epitopes: Population Coverage
	
	output$tblPeptides <- DT::renderDT ({
		# values$Active = "PP";
		# Reset 2nd & 3rd Tables:
		values$dfPopAlleles      = NULL;
		values$dfTotalPopulation = NULL;
		output$txtBtnDisplay = renderText(const$DisplayHLA$Warn);
		# Multiple Protein Sequences:
		nColTi = 8; # Col: Ti
		if(values$multSeq) {
			nmsCol = c("Peptide", "HLA", "Seq");
			nColTi = 9;
		} else {
			nmsCol = c("Peptide", "HLA");
		}
		dfPP = values$dfFltData[nmsCol];
		x = freq.all(dfPP, freq.hla(), type = values$typeHLA);
		values$dfPopCoverPP = x;
		#
		nmsCol = names.hla(values$typeHLA);
		DT::datatable(x, filter = 'top',
			options = option.regex(options$reg.PP,
				varia = list(order = list(nColTi, "desc")))) |>
		formatRound(c(nmsCol,'Tn','Ti'), 4);
	})
	
	### Overall Population: Cover & Alleles
	# Selected Rows:
	observeEvent(input$btnCovHLA, {
		ids = input$tblPeptides_rows_selected;
		print(ids); # DEBUG
		# NO Epitopes selected
		if(is.null(values$dfPopCoverPP) || length(ids) == 0) {
			values$dfPopAlleles  = NULL;
			values$dfTotalPopulation = NULL;
			values$fltHLAEpiSel  = NULL;
			output$txtBtnDisplay = renderText(const$DisplayHLA$Warn);
			return();
		}
		output$txtBtnDisplay = renderText(const$DisplayHLA$More);
		# Population Cover: Selected Epitopes
		pp = values$dfPopCoverPP[ids, "Peptide"];
		pp = unique(pp); # Avoid duplicates
		x  = values$dfFltData[c("HLA", "Peptide")];
		x  = x[x$Peptide %in% pp, ];
		x  = unique(x); # Avoid duplicates
		# HLA-Alleles covered:
		values$dfPopAlleles = x;
		# Population: Overall Coverage
		hla = sort(unique(x$HLA));
		xT  = freq.populationTotal(hla, options$HLA,
			type = values$typeHLA, digits = 3);
		values$dfTotalPopulation = xT;
		values$fltHLAEpiSel = hla;
	})
	
	output$tblAllelesPP = DT::renderDT(
			# old variant: dom = "lrtip"
			DT::datatable(values$dfPopAlleles, filter = 'top',
				options = option.regex(options$reg.PP,
					varia = list(dom = "tip", order = list(1, "asc"))))
	);
	output$tblTotalPopulation = DT::renderDT(
			DT::datatable(values$dfTotalPopulation, rownames = FALSE,
				options = option.regex(options$reg.PP,
					varia = list(dom = "t")))
	);
	
	# Pop Cover: Remaining HLA
	observeEvent(input$btnRemainingEpi, {
		if(is.null(values$dfFltData) ||
			nrow(values$dfFltData) == 0) return();
		values$pageTblPP = NULL;
		isHLA = values$dfFltData$HLA %in% values$fltHLAEpiSel;
		print(values$fltHLAEpiSel); # Debug
		cols  = if(values$multSeq) "Seq" else character(0);
		cols  = c("Peptide", "HLA", cols);
		dfHLA = values$dfFltData[! isHLA, cols, drop = FALSE];
		if(nrow(dfHLA) == 0) {
			print("NO more alleles!");
			output$tblRemainingEpi = NULL;
			values$dfRemainingEpi  = NULL;
			return();
		}
		# TODO: a lot of quasi-duplicated code;
		nColTi = 8; # Col: Ti
		# Multiple Protein Sequences:
		if(values$multSeq) {
			nColTi = 9;
		}
		# Filter:
		# https://datatables.net/reference/option/dom
		flt = getFilter.tblPopCovEpi();
		flt = list(order = list(nColTi, "desc"),
			searchCols = flt, dom = "tip");
		values$optRemainingEpi = option.regex(options$reg.PP, varia = flt);
		# Data:
		typeHLA = values$typeHLA;
		freqHLA = freq.population(dfHLA[c("Peptide", "HLA")],
					options$HLA, type = typeHLA);
		nmsCol  = cols;
		dfTmp   = dfHLA[nmsCol];
		x = freq.all(dfTmp, freqHLA, type = typeHLA);
		x = x[x$Ti > 0, , drop = FALSE]; # Exclude: HLA w. freq = 0;
		# Total Coverage:
		ids  = match(x$Peptide, values$dfPopCoverPP$Peptide);
		x$Tn = values$dfPopCoverPP$Ti[ids];
		idNm = match("Tn", names(x));
		names(x)[idNm] = "Ta"; # All Alleles (for Redundancy)
		# Decreasing Order:
		ids = order(x$Ti, x$Ta, decreasing = TRUE);
		x   = x[ids, ];
		values$dfRemainingEpi = x;
	})
	output$tblRemainingEpi = DT::renderDT({
		if(is.null(values$dfRemainingEpi)) return(NULL);
		nmsCol = names.hla(values$typeHLA);
		DT::datatable(values$dfRemainingEpi,
			filter = "top",
			options = values$optRemainingEpi) |>
		formatRound(c(nmsCol,'Ta','Ti'), 3);
	})
	
	getFilter.tblPopCovEpi = function() {
		flt   = input$tblPeptides_search_columns;
		nms   = if(values$multSeq) c("Len", "Seq") else "Len";
		idCol = match(nms, names(values$dfPopCoverPP));
		flt   = flt[idCol];
		if(all(flt == "")) {
			return(NULL);
		}
		nc = ncol(values$dfPopCoverPP) + 1; # Row id
		idCol  = idCol + 1;
		fltOld = lapply(flt, function(x) if(x == "") NULL else list(search = x));
		fltNew = rep(list(NULL), nc);
		for(id in seq_along(idCol)) fltNew[idCol[id]] = fltOld[id];
		return(fltNew);
	}
	
	output$txtPPTblPage = renderText({
		pg = values$pageTblPP;
		if(is.null(pg) || checkNoRemainingEpitopes())
			return("Remaining epitopes: not yet selected.");
		txt = paste0(pg$Epi, ": ", pg$Page, collapse = "; ");
		txt = paste0("Epitope on page: ", txt);
		return(txt);
	})
	
	# Remaining Epitope: Find Page
	checkNoRemainingEpitopes = function() {
		idR = input$tblRemainingEpi_rows_selected;
		if(is.null(idR) || length(idR) == 0)  {
			print("Reset page!");
			# values$pageTblPP = NULL;
			return(TRUE);
		}
		return(FALSE);
	}
	observeEvent(input$tblRemainingEpi_rows_selected, {
		if(checkNoRemainingEpitopes()) {
			values$pageTblPP = NULL;
			return();
		}
		# Selected Epitopes:
		idR = input$tblRemainingEpi_rows_selected;
		epi = values$dfRemainingEpi$Peptide[idR];
		# Main Table
		ide = input$tblPeptides_rows_all;
		tbl = values$dfPopCoverPP$Peptide[ide];
		idR = match(epi, tbl);
		pg  = (idR - 1) %/% 10 + 1; # TODO: Items per page;
		values$pageTblPP = list(Epi = epi, Page = pg);
	})
	
	# Remaining Epitopes: GoTo Page
	observeEvent(input$btnGoToPage, {
		pp = values$pageTblPP;
		if(is.null(pp) || length(pp$Page) == 0) return();
		# Debug:
		print(paste0("GoTo page: ", pp$Page[1]));
		proxy = dataTableProxy('tblPeptides');
		selectPage(proxy, pp$Page[1]);
	})
	
	### Save Data
	# Filtered Data:
	output$downloadData <- downloadHandler(
		filename = function() {
			paste("PP.Filtered", ".csv", sep = "");
		},
		content = function(file) {
			# TODO: check first if NULL;
			# All rows: https://rstudio.github.io/DT/shiny.html
			x = filter.byTable0();
			if(is.null(x)) return(NULL);
			write.csv(x, file, row.names = FALSE);
		}
	)
	
	# Download Epitopes:
	output$btnDownloadPP <- downloadHandler(
		filename = function() {
			paste("PP.Freq", ".csv", sep = "");
		},
		content = function(file) {
			# TODO: check first if NULL;
			# All rows: https://rstudio.github.io/DT/shiny.html
			ids = input$tblPeptides_rows_all;
			if(is.null(ids)) return(NULL);
			x = values$dfPopCoverPP[ids, ];
			write.csv(x, file, row.names = FALSE);
		}
	)
	
	### Print Selection
	observeEvent(input$printPPSel, {
		# TODO: check first if NULL;
		ids = input$tblPeptides_rows_selected;
		if(is.null(ids) || length(ids) == 0) return(NULL);
		x = values$dfPopCoverPP$Peptide[ids];
		print(x);
		showModal(modalDialog(
			title = "Selection",
			paste0(x, collapse = ", "),
			easyClose = TRUE, footer = NULL
		));
	})
	
	### Search Epitopes in Peptide
	# [Sub-Sequences]
	
	output$txtBtnSearchSubSeq = renderText({
		if(values$warnSubSeq) {
			txt = const$SearchNoSubSeq;
		} else {
			txt = const$SearchSubSeq;
		}
	});
	
	observeEvent(input$btnSearchSubSeq, {
		txt = input$inSubSeq;
		dat = if(options$allEpi.SubSeq) values$fullData else values$dfFltData;
		ids = find.subseq(txt, data = dat$Peptide, len = c(9, 11));
		if(is.null(ids) || length(ids) == 0) {
			output$tblSummarySubSeq = NULL;
			output$tblSubSeq  = NULL;
			values$warnSubSeq = TRUE;
			return();
		}
		values$warnSubSeq = FALSE;
		tbl = dat[ids,];
		tblSum = aggregate.subSeq(tbl);
		output$tblSummarySubSeq = renderDT(
			DT::datatable(tblSum, filter = 'top',
				options = option.regex(options$reg.PP)));
		output$tblSubSeq = renderDT(
			DT::datatable(tbl, filter = 'top',
				options = option.regex(options$reg.PP)));
	})
	
	
	### Summarize Epitopes
	
	output$txtBtnEpiSummary = renderText({
		if(values$warnEpiSummary) {
			txt = const$NoEpiSummary;
		} else {
			txt = const$EpiSummary;
		}
	});
	
	observeEvent(input$btnEpiSummary, {
		if(is.null(values$fullData)) return();
		txt = input$inEpiSummary;
		txt = strsplit(txt, "[, \n\t]+");
		txt = sort(txt[[1]]);
		# Process the entire Data set:
		dat = if(options$allEpi.EpiSummary) values$fullData else values$dfFltData;
		isRow = dat$Peptide %in% txt;
		dat = dat[isRow, ];
		if(nrow(dat) == 0) {
			output$tblEpiSummary  = NULL;
			values$warnEpiSummary = TRUE;
			return();
		}
		values$warnEpiSummary = FALSE;
		# TODO: summary
		output$tblEpiSummary = renderDT(
			DT::datatable(dat, filter = 'top',
				options = option.regex(options$reg.PP)));
	})
	
	### Protein Graph
	
	output$imgProtein <- renderPlot({
		nS = unique(values$dfFltData$start);
		plot.new();
		plot.window(xlim = range(nS), ylim = c(0, 2));
		for(npos in nS) {
			# lines(c(npos, npos), c(1,2), col = options$col.Pr, lwd = options$lwd.Pr);
			polygon(npos + c(-0.5, 0.5, 0.5, -0.5), c(1,1,2,2),
				col = options$col.Pr, border = options$border.Pr);
		}
	})
	
	### Help
	output$txtHelp <- renderUI({
		getHelp();
	})
}
