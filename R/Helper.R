

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


#################
### HLA Tools ###

trim.hla = function(x) {
	sub("(?i)^HLA-", "", x);
}

# Expand NULL to zero;
check.hla.df = function(x) {
	if(is.null(x$A)) x$A = 0;
	if(is.null(x$B)) x$B = 0;
	if(is.null(x$C)) x$C = 0;
	return(x);
}

### Tbl HLA Alleles
merge.hla = function(x, f, digits = 6) {
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

### Population Coverage
# x   = Set of Peptides w HLA-Alleles;
# hla = Frequency of HLA-Alleles;
freq.all = function(x, hla, seqPP = NULL, digits = 6) {
	if(is.null(x)) return(NULL);
	# Count(Alleles HLA)
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
	yAB  = y$A + y$B;
	y$Tn = yAB + y$C; # Total sum (naive total);
	y$Ti = round(y$Tn - yAB*y$C + y$A*y$B*(y$C - 1), digits);
	y$Tn = round(y$Tn, digits);
	return(y);
}

freq.populationTotal = function(x, f, do.totals = TRUE, digits = 6) {
	x = if(inherits(x, "data.frame")) x["HLA", drop = FALSE] else data.frame(HLA = x);
	x = merge(x, f, by = "HLA", all.x = TRUE);
	isMissing = is.na(x$Freq);
	x$Type[isMissing] = substr(x$HLA[isMissing], 1, 1);
	x$Freq[isMissing] = 0;
	tf  = tapply(x$Freq, x$Type, sum, na.rm = TRUE);
	nms = names(tf);
	tf = matrix(tf, nrow = 1);
	tf = round(tf, digits);
	tf = data.frame(tf); names(tf) = nms;
	tf = check.hla.df(tf);
	tf = tf[c("A", "B", "C")];
	#
	asZ = function(x) {
		if(is.na(x)) 0 else x;
	}
	tf$A = asZ(tf$A); tf$B = asZ(tf$B); tf$C = asZ(tf$C);
	if(do.totals) {
		tfAB  = tf$A + tf$B;
		Total = tfAB + tf$C;
		tf$Ti = round(Total - tfAB*tf$C - tf$A*tf$B*(1 - tf$C), digits);
		tf = cbind("Total" = "Population", tf); # "HLA" = c("Total")
	}
	return(tf);
}

####################

### Sub-Seq

#@@ Search for Epitopes in a given peptide
# len = range of sub-seq length;
find.subseq = function(x, data, len = c(9, 11), sort = TRUE, verbose = TRUE) {
	if(length(len) > 2) {
		if(verbose) cat("Only specified lengths will be used!\n");
	} else if(length(len) == 2) {
		len = sort(len);
		len = seq(len[1], len[2]);
	}
	if(all(nchar(data) < min(len))) return(NULL);
	# Generate all Sub-Seq:
	xSeq = subSeq(x, len, sort=sort);
	# Find matches:
	id = which(data %in% xSeq);
	return(id);
}

subSeq = function(x, len, sort = TRUE, unique = TRUE) {
	LEN  = nchar(x);
	sSeq = lapply(len, function(len) {
		if(LEN < len) return();
		sapply(seq(LEN - len + 1), function(nS) {
			substr(x, nS, nS + len - 1);
		})
	});
	sSeq = unlist(sSeq);
	if(unique) sSeq = unique(sSeq);
	if(sort)   sSeq = sort(sSeq);
	return(sSeq);
}

aggregate.subSeq = function(x, probs = 0.25, digits = 3, ...) {
	tblSum = tapply(x$Rank, x$Peptide, simplify = FALSE,
		function(x) {
			data.frame("", 0, min(x),
				round(quantile(x, probs=probs, na.rm = TRUE), digits = digits),
				length(x));
	});
	tblSum = do.call(rbind, tblSum);
	tblSum[1] = rownames(tblSum);
	tblSum[2] = nchar(tblSum[1]);
	rownames(tblSum) = NULL;
	nms = paste0("Rank", probs);
	nms = c("Peptide", "Len", "Rank", nms, "Count");
	names(tblSum) = nms;
	return(tblSum);
}
