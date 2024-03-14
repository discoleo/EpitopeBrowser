
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

freq.all = function(x, hla, digits = 5) {
	if(is.null(x)) return(NULL);
	x = data.frame(table(x));
	names(x)[1] = "Peptide";
	x$Peptide = as.character(x$Peptide);
	x$Len = nchar(x$Peptide);
	x = merge(x, hla, by = "Peptide");
	x$Total = round(x$A + x$B + x$C, digits);
	x$Ti = round(x$Total - (x$A + x$B)*x$C + x$A*x$B*(x$C - 1), digits);
	return(x);
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
