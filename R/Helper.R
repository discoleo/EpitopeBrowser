

### Basic Functions

### Read Epitopes
# sep = separator used in the csv file;
read.epi = function(file, hla.strip = TRUE, sep = ",") {
	x = read.csv(file, header = TRUE, sep = sep);
	if(is.numeric(x[,1])) {
		# IEDB: new format / unedited;
		hasMultiSeq = length(unique(x[,1])) > 1;
		id = if(hasMultiSeq) c(6,2,3,4,5, 1, 7,11,8) else c(6,2,3,4,5,7,11,8);
		x  = x[, id];
		if(hasMultiSeq) {
			idNms = c(1,2,5, 6,7,8,9);
			nms = c("HLA", "Peptide", "Len", "Seq", "ID", "Score","Rank");
		} else {
			idNms = c(1,2,5, 6,7,8);
			nms = c("HLA", "Peptide", "Len", "ID", "Score","Rank");
		}
		names(x)[idNms] = nms;
	} else if(is.numeric(x[,2])) {
		# Old / Modified csv-Files;
		names(x)[c(1,2,6,8)] = c("HLA", "Seq", "Peptide", "Rank");
	} else {
		# Saved Data: nothing;
	}
	if(hla.strip) x$HLA = trim.hla(x$HLA);
	return(x);
}


### HLA Tools

trim.hla = function(x) {
	sub("(?i)^HLA-", "", x);
}

check.hla.df = function(x) {
	if(is.null(x$A)) x$A = 0;
	if(is.null(x$B)) x$B = 0;
	if(is.null(x$C)) x$C = 0;
	return(x);
}

merge.hla = function(x, f, digits = 2) {
	if(nrow(x) == 0) warning("No data!");
	idHLA = match(c("HLA", "Freq"), names(f));
	f = f[, idHLA];
	f[, 2] = round(f[, 2] * 100, digits = digits);
	names(f)[2] = "Population (%)"
	x = merge(x, f, by = "HLA", all.x = TRUE);
	return(x);
}

freq.population = function(x, f) {
	if(nrow(x) == 0) warning("No data!");
	x  = merge(x[c("HLA", "Peptide")], f, by = "HLA", all.x = TRUE);
	isMissing = is.na(x$Freq);
	x$Type[isMissing] = substr(x$HLA[isMissing], 1, 1);
	x$Freq[isMissing] = 0;
	tf = tapply(x$Freq, x[c("Peptide", "Type")], sum, na.rm = TRUE);
	pp = rownames(tf);
	tf = data.frame(tf);
	tf$Peptide = pp;
	tf = check.hla.df(tf);
	#
	asZ = function(x) {
		isNA = is.na(x);
		x[isNA] = 0;
		return(x);
	}
	tf$A = asZ(tf$A); tf$B = asZ(tf$B); tf$C = asZ(tf$C);
	return(tf);
}

freq.all = function(x, hla, seqPP = NULL, digits = 5) {
	if(is.null(x)) return(NULL);
	y = data.frame(table(x));
	names(y)[1] = "Peptide";
	y$Peptide = as.character(y$Peptide);
	y$Len = nchar(y$Peptide);
	y = merge(y, hla, by = "Peptide");
	if(! is.null(seqPP)) {
		if(length(x) != length(seqPP)) {
			warning("Epitopes: Lengths do NOT match!");
		}
		id = match(y$Peptide, x);
		y$Seq = seqPP[id];
	}
	# Population Coverage:
	y$Total = round(y$A + y$B + y$C, digits);
	y$Ti = round(y$Total - (y$A + y$B)*y$C + y$A*y$B*(y$C - 1), digits);
	return(y);
}

freq.populationTotal = function(x, f, digits = 3) {
	x = if(inherits(x, "data.frame")) x["HLA", drop = FALSE] else data.frame(HLA = x);
	x = merge(x, f, by = "HLA", all.x = TRUE);
	isMissing = is.na(x$Freq);
	x$Type[isMissing] = substr(x$HLA[isMissing], 1, 1);
	x$Freq[isMissing] = 0;
	tf  = tapply(x$Freq, x$Type, sum, na.rm = TRUE);
	nms = names(tf);
	tf  = matrix(tf, nrow = 1);
	tf = round(tf, digits);
	tf = data.frame(tf); names(tf) = nms;
	tf = check.hla.df(tf);
	tf = tf[c("A", "B", "C")];
	#
	asZ = function(x) {
		if(is.na(x)) 0 else x;
	}
	tf$A = asZ(tf$A); tf$B = asZ(tf$B); tf$C = asZ(tf$C);
	return(tf);
}
