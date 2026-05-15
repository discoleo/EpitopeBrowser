

#################
### HLA Tools ###

### Trim HLA
trim.hla = function(x) {
	sub("(?i)^HLA-", "", x);
}

### Remove HLA_DX A-Chain
trim.hla2A = function(x, trimmed = TRUE) {
	if(trimmed) {
		sub("(?i)^D[PQR]A[1-5][0-9*\\:]++/", "", x, perl = TRUE);
	} else {
		sub("(?i)^HLA-D[PQR]A[1-5][0-9*\\:]++/", "", x, perl = TRUE);
	}
}


### HLA Class
classHLA = function(x, stripHLA = TRUE) {
	if(stripHLA) {
		hla2 = "(?i)^D";
	} else {
		hla2 = "(?i)^HLA-D";
	}
	isHLA2 = grepl(hla2, x$HLA);
	if(all(isHLA2)) {
		return(2);
	} else if(any(isHLA2)) {
		return(0); # TODO
	}
	return(1);
}

### Collapse Alleles
#' @export
collapse.hla = function(x, hla.add, sep = ", ") {
	if(hla.add) {
		x$HLA = paste0("HLA-", x$HLA);
	}
	x = tapply(x$HLA, x$Peptide, paste0, collapse = sep);
	x = cbind(names(x), x);
	rownames(x) = NULL;
	return(x);
}

### Std Names
names.hla = function(type = 1) {
	if(type == 1) c("A", "B", "C")
		else c("DP", "DQ", "DR", "DR3");
}
names.hlaFreq = function(type = 1) {
	c(names.hla(type=type), "Tn", "Ti");
}


### Filter HLA Alleles
# dfHLA = Set of HLA Alleles w. Frequencies;
filter.HLA = function(x, fltAllele, dfHLA, hla.strip = TRUE) {
	if(fltAllele == "All") return(x);
	# Data with HLA:
	hla = if(hla.strip) x$HLA else trim.hla(x$HLA);
	if(fltAllele == ">= 3%") {
		frequentHLA = dfHLA$HLA[dfHLA$Freq >= 0.03];
		isHLA = hla %in% frequentHLA;
	} else if(fltAllele == ">= 1%") {
		frequentHLA = dfHLA$HLA[dfHLA$Freq >= 0.01];
		isHLA = hla %in% frequentHLA;
	} else if(fltAllele == "< 1%") {
		frequentHLA = dfHLA$HLA[dfHLA$Freq < 0.01];
		isHLA = hla %in% frequentHLA;
	} else {
		isHLA = hla %in% dfHLA$HLA;
		if(fltAllele == "Uncommon") isHLA = ! isHLA;
	}
	# print("HLA"); print(head(dfHLA))
	x = x[isHLA, ];
	return(x);
}

