

### Basic Functions

### Strings

### Text => Values of Range
splitRange = function(x) {
	txt = strsplit(x, "[-, /]+")[[1]];
	val = round(as.numeric(txt));
	val = val[! is.na(val)];
	return(val);
}

### Read Epitopes
# sep = separator used in the csv file;
read.epi = function(file, hla.strip = TRUE, sep = ",") {
	x = read.csv(file, header = TRUE, sep = sep);
	if(is.numeric(x[,1])) {
		# IEDB: new format / unedited;
		hasMultiSeq = length(unique(x[,1])) > 1;
		# MHC-2 vs MGC-1:
		idMHC2 = which(grepl("netmhciipan_el.score", names(x)));
		isMHC2 = length(idMHC2) > 0;
		#
		id = if(isMHC2) idMHC2 else 11;
		id = if(hasMultiSeq) c(6,2,3,4,5, 1, 7,id,8) else c(6,2,3,4,5,7,id,8);
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

trim.hla2A = function(x) {
	sub("(?i)^D[PQR]A[1-5][0-9*\\:]++/", "", x, perl = TRUE);
}

### Std Names
names.hla = function(type = 1) {
	if(type == 1) c("A", "B", "C")
		else c("DP", "DQ", "DR");
}

### FIlter HLA Alleles
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


### Checks & Corrections

# Expand NULL to zero;
check.hla.df = function(x, type = 1) {
	if(type == 1) {
		if(is.null(x$A)) x$A = 0;
		if(is.null(x$B)) x$B = 0;
		if(is.null(x$C)) x$C = 0;
	} else {
		if(is.null(x$DP)) x$DP = 0;
		if(is.null(x$DQ)) x$DQ = 0;
		if(is.null(x$DR)) x$DR = 0;
	}
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

### Population Coverage: each Epitope
# f = HLA Frequency;
freq.population = function(x, f, type = 1) {
	if(nrow(x) == 0) warning("No data!");
	x = unique(x[c("HLA", "Peptide")]);
	x = merge(x, f, by = "HLA", all.x = TRUE);
	# Missing Alleles:
	isMissing = is.na(x$Freq);
	LEN.HLA = if(type == 1) 1 else 2; # TODO: robust; full HLA-2?
	x$Type[isMissing] = substr(x$HLA[isMissing], 1, LEN.HLA);
	x$Freq[isMissing] = 0;
	tf = tapply(x$Freq, x[c("Peptide", "Type")], sum, na.rm = TRUE);
	pp = rownames(tf);
	tf = data.frame(tf);
	tf$Peptide = pp;
	tf = check.hla.df(tf, type=type);
	asZ = function(x) {
		isNA = is.na(x);
		x[isNA] = 0;
		return(x);
	}
	if(type == 1) {
		tf$A = asZ(tf$A); tf$B = asZ(tf$B); tf$C = asZ(tf$C);
	} else {
		tf$DP = asZ(tf$DP); tf$DQ = asZ(tf$DQ);
		tf$DR = asZ(tf$DR);
	}
	return(tf);
}

### Population Coverage
# x   = Vector of Peptides (PP repeated for each HLA-Allele);
# hla = List of Epitopes w. Population Cover;
# seqPP = Protein Seq;
freq.all = function(x, hla, seqPP = NULL, type = 1, digits = 6) {
	if(is.null(x)) return(NULL);
	# Count(Alleles HLA)
	isDF = inherits(x, "data.frame");
	if(isDF) {
		noSeq = is.null(x$Seq);
		seqPP = if(noSeq) NULL else unique(x[c("Peptide", "Seq")]);
		x = unique(x[c("Peptide", "HLA")]);
		x = x$Peptide;
	}
	y = data.frame(table(x));
	names(y)[1] = "Peptide";
	y$Peptide = as.character(y$Peptide);
	y$Len = nchar(y$Peptide);
	y = merge(y, hla, by = "Peptide");
	# Add Seq-ID:
	if(! is.null(seqPP)) {
		if(isDF) {
			y = merge(y, seqPP, by = "Peptide");
		} else {
			if(length(x) != length(seqPP)) {
				warning("Epitopes: Lengths do NOT match!");
			}
			id = match(y$Peptide, x);
			y$Seq = seqPP[id];
		}
	}
	# Population Coverage:
	if(type == 1) {
		yA = y$A; yB = y$B; yC = y$C;
	} else {
		yA = y$DP; yB = y$DQ; yC = y$DR;
	}
	yAB  = yA  + yB;
	y$Tn = yAB + yC; # Total sum (naive total);
	y$Ti = round(y$Tn - yAB*yC + yA*yB*(yC - 1), digits);
	y$Tn = round(y$Tn, digits);
	return(y);
}

### Total Coverage
# (for ALL Epitopes)
# x = HLA Alleles;
# f = Freq of HLA Alleles;
# type = HLA class 1 vs class 2;
freq.populationTotal = function(x, f, type = 1, do.totals = TRUE, digits = 6) {
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
	tf = check.hla.df(tf, type=type);
	tf = tf[names.hla(type)]; # Class-1 vs Class-2
	#
	asZ = function(x) {
		if(is.na(x)) 0 else x;
	}
	if(type == 1) {
		tfA = asZ(tf$A); tf$A = tfA;
		tfB = asZ(tf$B); tf$B = tfB;
		tfC = asZ(tf$C); tf$C = tfC;
	} else {
		tfA = asZ(tf$DP); tf$DP = tfA;
		tfB = asZ(tf$DQ); tf$DQ = tfB;
		tfC = asZ(tf$DR); tf$DR = tfC;
	}
	if(do.totals) {
		tfAB  = tfA  + tfB;
		Total = tfAB + tfC;
		tf$Ti = round(Total - tfAB*tfC - tfA*tfB*(1 - tfC), digits);
		tf = cbind("Total" = "Population", tf); # "HLA" = c("Total")
	}
	return(tf);
}

####################

### HLA Regions

merge.RegionsHLA = function() {
	nms = c("HLA", "Freq");
	sfx = paste0("Freq.", c("De", "It", "Hu"));
	x = hla("De");
	y = hla("Italy")[, nms];
	names(x)[3] = sfx[1];
	names(y)[2] = sfx[2];
	x = merge(x, y, by = "HLA", all = TRUE);
	# Level 1:
	x$HLA.L1 = sub("\\:[0-9]++$", "", x$HLA, perl = TRUE);
	y = hla("Hungary")[, nms];
	y = y[! is.na(y$HLA), ];
	# Exclude Level 2:
	y = y[! grepl("\\:[0-9]++$", y$HLA, perl = TRUE), ];
	names(y) = c("HLA.L1", sfx[3]);
	x = merge(x, y, by = "HLA.L1", all = TRUE);
	print(str(x))
	return(x);
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
