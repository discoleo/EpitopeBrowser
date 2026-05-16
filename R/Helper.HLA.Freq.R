
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

### Tbl HLA Alleles
# y = frequency of individual HLA alleles;
#' @export
merge.hla = function(x, y, digits = 6) {
	if(nrow(x) == 0) warning("No data!");
	idHLA = match(c("HLA", "Freq"), names(y));
	f = y[, idHLA];
	# Round:
	if(! is.null(digits)) {
		f[, 2] = round(f[, 2] * 100, digits = digits);
		names(f)[2] = "Population (%)";
	} else names(f)[2] = "Population (pr)";
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
	# TODO: robust; full HLA-2?
	LEN.HLA = if(type == 1) 1 else 2;
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
		tf$DR3 = asZ(tf$DR3);
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
	y = popCover(y, type=type, digits=digits);
	return(y);
}

popCover = function(x, type, digits = NULL) {
	y = x; yDR3 = 0;
	if(type == 1) {
		yA = y$A; yB = y$B; yC = y$C;
	} else {
		yA = y$DP; yB = y$DQ; yC = y$DR;
		yDR3 = y$DR3;
	}
	yAB  = yA  + yB;
	y$Tn = yAB + yC + yDR3; # Total sum (naive total);
	# y$Ti = y$Tn - yAB*yC + yA*yB*(yC - 1); # [old]
	y$Ti = 1 - (1-yA)*(1-yB)*(1-yC)*(1-yDR3);
	if(! is.null(digits)) {
		y$Ti = round(y$Ti, digits);
		y$Tn = round(y$Tn, digits);
	}
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
	tfDR3 = 0;
	if(type == 1) {
		tfA = asZ(tf$A); tf$A = tfA;
		tfB = asZ(tf$B); tf$B = tfB;
		tfC = asZ(tf$C); tf$C = tfC;
	} else {
		tfA = asZ(tf$DP); tf$DP = tfA;
		tfB = asZ(tf$DQ); tf$DQ = tfB;
		tfC = asZ(tf$DR); tf$DR = tfC;
		tfDR3 = asZ(tf$DR3); tf$DR3 = tfDR3;
	}
	if(do.totals) {
		# [old]
		# tfAB  = tfA  + tfB;
		# Total = tfAB + tfC;
		# tf$Ti = round(Total - tfAB*tfC - tfA*tfB*(1 - tfC), digits);
		tf$Ti = round(1 - (1-tfA)*(1-tfB)*(1-tfC)*(1-tfDR3), digits);
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
freq.loc2any = function(p1, p2) {
	n1 = length(p1); n2 = length(p2);
	p  = 1 - p1[n1]^2 * p2[n2]^2;
	return(p);
}
freq.loc2any.old = function(p1, p2) {
	# Step by step computation;
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
		if(is.null(x$DR3)) x$DR3 = 0;
	}
	return(x);
}
