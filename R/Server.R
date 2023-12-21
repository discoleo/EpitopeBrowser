

getServer = function(x) {
	shiny::shinyServer(server.app)
}

server.app = function(input, output, session) {
	# Dynamic variable
    values = reactiveValues();
	
	options = list(
		hla.strip = TRUE, # HLA-A... => A...;
		sep = ",",       # csv Separator
		HLA = as.data.frame.hla(),
		reg.Data = TRUE,  # Regex for Data-Table
		reg.PP   = TRUE   # Regex for Epitopes-Table
	);
	
	freq.hla = function() {
		x = values$fData[c("Peptide", "HLA")];
		# isDuplic = duplicated(x);
		# x = x[ ! isDuplic, ];
		freq.population(x, options$HLA);
	}
	
	trim = function(x) {
		sub("(?i)^HLA-", "", x);
	}
	filter.df = function(x, lim.rank, fltAllele) {
		x = x[x$Rank <= lim.rank, ];
		x = filter.HLA(x, fltAllele);
		return(x);
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
	#
	option.regex = function(x, caseInsens = TRUE, varia = NULL) {
		opt = list(search = list(regex = x, caseInsensitive = caseInsens));
		if( ! is.null(varia)) opt = c(opt, varia);
		return(opt);
	}
    # File Input
    observe({
        file1 <- input$file;
        if (is.null(file1)) 
            return(NULL);
        
		x = read.csv(file1$datapath, header = TRUE, sep = options$sep);
		names(x)[c(1,2,6,8)] = c("HLA", "Seq", "Peptide", "Rank");
		if(options$hla.strip) x$HLA = trim(x$HLA);
		#
        values$mydata = x;
    })
	
	### Options: Data
	observe({
		lim.rank = input$rank;
		fltAllele = input$fltAllele;
        x = values$mydata;
		values$fData = filter.df(x, lim.rank, fltAllele);
		print(nrow(values$fData))
	})
	observe({
		isReg = input$chkRegex;
		options$reg.Data = isReg;
	})
	
	### Tables
	output$tblData <- DT::renderDT ({
		DT::datatable(values$fData, filter = 'top',
			options = option.regex(options$reg.Data));
	})
	output$tblAlleles <- DT::renderDT ({
		x = values$fData$HLA;
		x = data.frame(table(x));
		names(x)[1] = "HLA";
		DT::datatable(x, filter = 'top');
	})
	output$tblPeptides <- DT::renderDT ({
		x = data.frame(table(values$fData$Peptide));
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
			return();
		}
		pp = values$pp[ids, "Peptide"];
		x  = values$fData[c("HLA", "Peptide")]
		x  = x[x$Peptide %in% pp, ];
		output$tblAllelesPP = DT::renderDT(
			DT::datatable(x, filter = 'top',
				options = option.regex(options$reg.PP,
					varia = list(dom = "lrtip", order = list(1, "asc"))))
		);
	})
}
