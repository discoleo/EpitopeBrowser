

getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		fltRank   = 0.55, # Filter: only <= limit;
		hla.strip = TRUE, # HLA-A... => A...;
		sep = ",",        # csv Separator
		HLA = as.data.frame.hla(),
		reg.Data  = TRUE,  # Regex for Data-Table
		reg.PP    = TRUE,  # Regex for Epitopes-Table
		highlight = TRUE,  # Highlight Search Term
		col.Pr    = "#FF0032A0",
		border.Pr = "#640000A0",
		lwd.Pr   = 1.5 # NOT used with Polygon
	);
	
	updateNumericInput(session, "rank", value = options$fltRank);
	
	const = list(Warn = list(
		DisplayHLA = "Select first some Epitopes in the table below.")
	);
	
	# Dynamic variable
	values = reactiveValues(
		Active    = "Not",  # Active Tab
		fullData  = NULL,   # initial Data
		dfGlData  = NULL,   # Globally filtered Data
		dfFltData = NULL,   # Data filtered in Table
		reg.Data  = options$reg.Data,
		fltRank   = NULL,   # is set automatically
		fltAllele = NULL,
		fltCols   = NULL, # TODO
		multSeq   = FALSE
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
	
	# Epitopes
	output$tblPeptides <- DT::renderDT ({
		# values$Active = "PP";
		output$tblAllelesPP = NULL; # Reset 2nd & 3rd Tables;
		output$tblTotalPopulation = NULL;
		output$txtBtnDisplay = renderText(const$Warn$DisplayHLA);
		# Multiple Protein Sequences:
		if(values$multSeq) {
			idSeq = values$dfFltData$Seq;
		} else idSeq = NULL;
		x = freq.all(values$dfFltData$Peptide, freq.hla(), seqPP = idSeq);
		values$pp = x;
		#
		DT::datatable(x, filter = 'top',
			options = option.regex(options$reg.PP,
				# nCol = 8: Ti;
				varia = list(order = list(8, "desc"))));
	})
	
	# Selected Rows:
	observeEvent(input$ppHLA, {
		ids = input$tblPeptides_rows_selected;
		print(ids); # DEBUG
		# NO Epitopes selected
		if(is.null(values$pp) || length(ids) == 0) {
			output$tblAllelesPP = NULL;
			output$tblTotalPopulation = NULL;
			output$txtBtnDisplay = renderText(const$Warn$DisplayHLA);
			return();
		}
		output$txtBtnDisplay = NULL;
		pp = values$pp[ids, "Peptide"];
		x  = values$dfFltData[c("HLA", "Peptide")]
		x  = x[x$Peptide %in% pp, ];
		output$tblAllelesPP = DT::renderDT(
			DT::datatable(x, filter = 'top',
				options = option.regex(options$reg.PP,
					varia = list(dom = "lrtip", order = list(1, "asc"))))
		);
		# Total Population:
		hla = unique(x$HLA);
		xT = freq.populationTotal(hla, options$HLA);
		Total = xT$A + xT$B + xT$C;
		xT$Ti = round(Total - (xT$A + xT$B)*xT$C - xT$A*xT$B + xT$A*xT$B*xT$C, 3);
		xT = cbind("Total" = "Population", xT);
		output$tblTotalPopulation = DT::renderDT(
			DT::datatable(xT, rownames = FALSE,
				options = option.regex(options$reg.PP,
					varia = list(dom = "t")))
		);
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
