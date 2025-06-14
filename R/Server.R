

getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		fltRank   = 0.55, # Default value for Rank-Filter;
		hla.strip = TRUE, # HLA-A... => A...;
		sep = ",",        # csv Separator
		HLA = as.data.frame.hla(),
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
	
	const = list(Warn = list(
		DisplayHLA     = "Select first some Epitopes in the table below."),
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
		multSeq   = FALSE,
		# Population Coverage:
		dfAllelesPP = NULL,
		dfTotalPopulation = NULL,
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
	reset.tab = function() {}
	
	freq.hla = function() {
		x = values$dfFltData[c("Peptide", "HLA")];
		freq.population(x, options$HLA);
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
		x = x[x$Rank <= lim.rank, ];
		x = filter.HLA(x, fltAllele);
		#
		values$dfGlData = x;
		cat("Rows: ", nrow(x), "\n");
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
	dataTable = function() ({
		# print("Rendering table!");
		flt = values$fltCols; # getFilter.tblData();
		if(! is.null(flt)) flt = list(searchCols = flt);
		DT::datatable(values$dfGlData, filter = 'top',
			options = option.regex(values$reg.Data, varia = flt));
	})
	
	output$tblData = DT::renderDT(dataTable())
	
	# HLA Alleles
	output$tblAlleles <- DT::renderDT ({
		# values$Active = "Alleles";
		x = values$dfFltData$HLA;
		x = data.frame(table(x));
		names(x)[1] = "HLA";
		x = merge.hla(x, options$HLA);
		DT::datatable(x, filter = 'top');
	})
	
	### Epitopes: Population Coverage
	
	output$tblPeptides <- DT::renderDT ({
		# values$Active = "PP";
		# Reset 2nd & 3rd Tables:
		values$dfAllelesPP = NULL;
		values$dfTotalPopulation = NULL;
		output$txtBtnDisplay = renderText(const$Warn$DisplayHLA);
		# Multiple Protein Sequences:
		nColTi = 8; # Col: Ti
		if(values$multSeq) {
			idSeq = values$dfFltData$Seq;
			nColTi = 9;
		} else idSeq = NULL;
		x = freq.all(values$dfFltData$Peptide, freq.hla(), seqPP = idSeq);
		values$pp = x;
		#
		DT::datatable(x, filter = 'top',
			options = option.regex(options$reg.PP,
				varia = list(order = list(nColTi, "desc"))));
	})
	
	# Selected Rows:
	observeEvent(input$btnCovHLA, {
		ids = input$tblPeptides_rows_selected;
		print(ids); # DEBUG
		# NO Epitopes selected
		if(is.null(values$pp) || length(ids) == 0) {
			values$dfAllelesPP = NULL;
			values$dfTotalPopulation = NULL;
			output$txtBtnDisplay = renderText(const$Warn$DisplayHLA);
			return();
		}
		output$txtBtnDisplay = NULL;
		pp = values$pp[ids, "Peptide"];
		x  = values$dfFltData[c("HLA", "Peptide")]
		x  = x[x$Peptide %in% pp, ];
		# HLA-Alleles covered:
		values$dfAllelesPP = x;
		# Population: Overall Coverage
		hla = unique(x$HLA);
		xT  = freq.populationTotal(hla, options$HLA, digits = 3);
		values$dfTotalPopulation = xT;
	})
	
	output$tblAllelesPP = DT::renderDT(
			DT::datatable(values$dfAllelesPP, filter = 'top',
				options = option.regex(options$reg.PP,
					varia = list(dom = "lrtip", order = list(1, "asc"))))
	);
	output$tblTotalPopulation = DT::renderDT(
			DT::datatable(values$dfTotalPopulation, rownames = FALSE,
				options = option.regex(options$reg.PP,
					varia = list(dom = "t")))
	);
	
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
	# Epitopes:
	output$downloadPP <- downloadHandler(
		filename = function() {
			paste("PP.Freq", ".csv", sep = "");
		},
		content = function(file) {
			# TODO: check first if NULL;
			# All rows: https://rstudio.github.io/DT/shiny.html
			ids = input$tblPeptides_rows_all;
			if(is.null(ids)) return(NULL);
			x = values$pp[ids, ];
			write.csv(x, file, row.names = FALSE);
		}
	)
	
	### Print Selection
	observeEvent(input$printPPSel, {
		# TODO: check first if NULL;
		ids = input$tblPeptides_rows_selected;
		if(is.null(ids) || length(ids) == 0) return(NULL);
		x = values$pp$Peptide[ids];
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
		HTML("<div>
			<p><b>Rank</b></p>
			IEDB recommends a cutoff value of <b>0.5</b> (for HLA-1 epitopes).
			A higher cutoff value (up to <b>1 - 1.5</b>) may still work well for epitopes derived from viral structural proteins. A stringent cutoff value (e.g. <b>0.2</b>) may be only rarely required.
			</div>"
		)
	})
}
