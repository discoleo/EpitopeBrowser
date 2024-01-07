

getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Global Options
	options = list(
		hla.strip = TRUE, # HLA-A... => A...;
		sep = ",",       # csv Separator
		HLA = as.data.frame.hla(),
		reg.Data = TRUE,  # Regex for Data-Table
		reg.PP   = TRUE   # Regex for Epitopes-Table
	);
	
	# Dynamic variable
	values = reactiveValues(
		Active    = "Not",  # Active Tab
		fullData  = NULL,   # initial Data
		fltGlData = NULL,   # Globally filtered Data
		fltData   = NULL,   # Data filtered in Table
		reg.Data  = options$reg.Data
	);
	
	# Reset Filters on Data table
	hasData = function() { ! is.null(values$fltGlData); }
	reset.tab = function() {
		if(values$Active != "Data" && hasData()) {
			values$Active  = "Data";
			values$fltData = values$fltGlData;
		}
	}
	
	observeEvent(input$menu.top, {
		if(values$Active == "Data") {
			print("Switched");
			filter.byTable();
			values$Active = "Other";
		}
	})
	
	freq.hla = function() {
		x = values$fltData[c("Peptide", "HLA")];
		# isDuplic = duplicated(x);
		# x = x[ ! isDuplic, ];
		freq.population(x, options$HLA);
	}
	
	trim = function(x) {
		sub("(?i)^HLA-", "", x);
	}
	filter.df = function() {
		lim.rank  = input$rank;
		fltAllele = input$fltAllele;
		x = values$fullData;
		x = x[x$Rank <= lim.rank, ];
		x = filter.HLA(x, fltAllele);
		#
		values$fltGlData = x;
		values$fltData   = x;
		cat("Rows: ", nrow(x), "\n");
		output$tblData = DT::renderDT(dataTable());
	}
	filter.HLA = function(x, fltAllele) {
		if(fltAllele == "all") return(x);
		#
		hla = if(options$hla.strip) x$HLA else trim(x$HLA);
		isHLA = hla %in% options$HLA$HLA;
		if(fltAllele == "uncommon") isHLA = ! isHLA;
		# print("HLA"); print(head(options$HLA))
		x = x[isHLA, ];
		return(x);
	}
	filter.byTable = function() {
		id = input$tblData_rows_all;
		values$fltData = values$fltGlData[id, ];
	}
	#
	option.regex = function(x, varia = NULL, caseInsens = TRUE) {
		opt = list(search = list(regex = x, caseInsensitive = caseInsens));
		if( ! is.null(varia)) opt = c(opt, varia);
		return(opt);
	}
	
	# File Input
	observe({
		file1 = input$file;
		if (is.null(file1)) 
			return(NULL);
		
		x = read.csv(file1$datapath, header = TRUE, sep = options$sep);
		names(x)[c(1,2,6,8)] = c("HLA", "Seq", "Peptide", "Rank");
		if(options$hla.strip) x$HLA = trim(x$HLA);
		#
		values$fullData = x;
	})
	
	### Options: Data
	observeEvent(values$fullData, {
		filter.df();
	})
	observeEvent(input$rank, {
		filter.df();
	})
	observeEvent(input$fltAllele, {
		filter.df();
	})
	observeEvent(input$chkRegex, {
		isReg = input$chkRegex;
		values$reg.Data = isReg;
		flt = input$tblData_search_columns;
		flt = c("", flt); # Row ID
		flt = lapply(flt, function(x) if(x == "") NULL else list(search = x));
		output$tblData <- DT::renderDT(dataTable(varia = list(searchCols = flt)));
	})
	
	# observeEvent(input$tblData_search_columns, {
	#	filter.byTable();
	# })
	
	### Tables
	
	# Data
	dataTable = function(varia = NULL) ({
		reset.tab();
		DT::datatable(values$fltGlData, filter = 'top',
			options = option.regex(values$reg.Data, varia=varia));
	})
	# probably NOT needed
	# output$tblData <- DT::renderDT ({
	#	reset.tab();
	#	DT::datatable(values$fltGlData, filter = 'top',
	#		options = option.regex(values$reg.Data));
	# })
	
	# HLA Alleles
	output$tblAlleles <- DT::renderDT ({
		# values$Active = "Alleles";
		x = values$fltData$HLA;
		x = data.frame(table(x));
		names(x)[1] = "HLA";
		DT::datatable(x, filter = 'top');
	})
	
	# Epitopes
	output$tblPeptides <- DT::renderDT ({
		# values$Active = "PP";
		output$tblAllelesPP = NULL; # Reset 2nd & 3rd Tables;
		output$tblTotalPopulation = NULL;
		x = data.frame(table(values$fltData$Peptide));
		names(x)[1] = "Peptide";
		x$Peptide = as.character(x$Peptide);
		x$Len = nchar(x$Peptide);
		x = merge(x, freq.hla(), by = "Peptide");
		x$Total = round(x$A + x$B + x$C, 5);
		x$Ti = round(x$Total - (x$A + x$B)*x$C + x$A*x$B*(x$C - 1), 5);
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
		print(ids);
		if(is.null(values$pp) || length(ids) == 0) {
			output$tblAllelesPP = NULL;
			output$tblTotalPopulation = NULL;
			return();
		}
		pp = values$pp[ids, "Peptide"];
		x  = values$fltData[c("HLA", "Peptide")]
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
}
