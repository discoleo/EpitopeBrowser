

### Basic Functions

### Strings

### Text => Values of Range
splitRange = function(x) {
	txt = strsplit(x, "[-, /]+")[[1]];
	val = round(as.numeric(txt));
	val = val[! is.na(val)];
	return(val);
}

table.length = function(x, hasMultiSeq, unique = TRUE) {
	colSeq = if(hasMultiSeq) "Seq" else NULL;
	x = x[, c("Peptide", colSeq, "start", "Len")];
	# Count each epitope only once;
	if(unique) {
		x = unique(x);
	}
	tbl = table(x$Len);
	nms = names(tbl);
	tbl = cbind("Count:", matrix(format(tbl), nrow = 1));
	tbl = as.data.frame(tbl);
	names(tbl) = c("Length:", nms);
	return(tbl);
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
names.hlaFreq = function(type = 1) {
	c(names.hla(type=type), "Tn", "Ti");
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


### Frequency: in 1 Step
# hla  = df with HLA-Frequencies;
# type = Type of HLA;
freq.hlaMerge = function(x, hla, multiSeq, type = 1, probs = c(0, 0.5, 1)) {
	colSeq   = if(multiSeq) "Seq" else NULL;
	datFreq  = freq.population(x, hla, type = type);
	datFreq  = freq.all(x[c("Peptide", "HLA", colSeq)], datFreq, type = type);
	datFreq$Len = NULL;
	datStats = summary.epi(x, multiSeq, probs = probs);
	datStats = merge(datStats, datFreq, by = c("Peptide", colSeq));
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

### Diploid Genome: Exact Frequency

# Note:
# - does NOT account for linkage disequilibrium;
# - computes coverage properly for Diploid genome;

### 2 Loci: 2 Alleles
# p1 = Frequency of specific allele at Locus 1;
# p2 = Frequency of (another) specific allele at Locus 2;
freq.loc2a2 = function(p1, p2) {
	p10 = 1 - p1; p20 = 1 - p2;
	# p1m = p1*p10; p2m = p2*p20;
	# p = p1^2 + p2^2 - p1^2*p2^2 + 4*p1m*p2m + 2*p1m*p20^2 + 2*p2m*p10^2;
	# Direct formula:
	p = 1 - p10^2 * p20^2;
	return(p);
}

### 2 Loci: Variable Alleles
# Last allele: considered 0;
# TODO: think more thoroughly the 0 case;
freq.loc2an = function(p1, p2) {
	n1 = length(p1);
	n2 = length(p2);
	g  = expand.grid(seq(p1), seq(p1), seq(p2), seq(p2));
	last = nrow(g);
	g = g[- last, ]; # Case: 0;
	tmp = apply(g, 1, \(id) {
		prod(c(p1[id[1:2]], p2[id[3:4]]));
	});
	sum(tmp);
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

names.quant = function(x) {
	names(x)[grepl("^Rq", names(x))];
}

summary.epi = function(x, isMultiSeq = NULL, probs = c(0, 0.5, 1)) {
	if(is.null(isMultiSeq)) {
		isMultiSeq = match("Seq", names(x));
		isMultiSeq = ! is.na(isMultiSeq);
	}
	# Unique Set:
	datUnq   = unique(x[c("Peptide", "Seq")]);
	datStats = tapply(x$Rank, x$Peptide, function(x) {
		as.data.frame(as.list(quantile(x, probs = probs)));
	});
	pp = names(datStats);
	datStats = do.call(rbind, datStats);
	names(datStats)  = paste0("Rq", format(probs));
	datStats$Peptide = pp;
	rownames(datStats) = NULL;
	#
	datUnq$Len = nchar(datUnq$Peptide);
	datStats   = merge(datUnq, datStats, by = "Peptide");
	return(datStats);
}
