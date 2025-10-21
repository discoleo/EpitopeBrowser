
### pI:
# Synthetic peptides used as vaccines should have a pI between 5 and 7.
# Peptides with pI below 5 or above 8 should probably be avoided.


### Theoretical pI:
# Ref:
# 1. R Package Peptides: Calculate Indices and Theoretical Physicochemical Properties
#    of Protein Sequences.
#    https://cran.r-project.org/web/packages/Peptides/index.html


# x = Vector of Peptide-Sequences;
# type  = set of pKa to use;
# split = aa-sequence needs to be split into individual AAs;
#' @export
pI = function(x, type = c("EMBOSS", "Bjellqvist", "Dawson",
		"Lehninger", "Murray"),
		split = TRUE, verbose = FALSE) {
	type = match.arg(type);
	pKa  = aa.pK(type);
	# Split into sequence of AAs;
	if(split) x = strsplit(x, "", fixed = TRUE);
	#
	aaCharge = function(pH) {
		div = 1.0 + 10^(pKa$charges * (pH - pKa$pK));
		charge = pKa$charges / div;
	}
	absoluteCharge = function(x, seqAA) {
		# x = pH;
		charge = aaCharge(x);
		charge = sum(charge * seqAA);
		charge = abs(charge); # for Optimisation!
		return(charge);
	}
	# AAs:
	LEN = length(pKa$AA);
	AAs = pKa$AA[ - c(LEN-1, LEN)];
	#
	pI = lapply(x, function(x) {
		fAA = table(x);
		ids = match(AAs, names(fAA));
		fAA = ifelse(is.na(ids), 0, fAA[ids]);
		fAA = c(fAA, 1, 1); # add: ('c', 'n')
		#
		pI = optimize(f = absoluteCharge, interval = c(0,14), seqAA = fAA);
		if(verbose) {
			pI = data.frame(pI = pI[[1]], charge = pI[[2]]);
			return(pI);
		}
		pI = unlist(pI)[[1]];
	});
	#
	if(verbose) {
		pI = do.call(rbind, pI);
	} else pI = unlist(pI);
	return(pI);
}


### Electric Charge at given pH

#' @export
charge = function(pH, x, type = c("EMBOSS", "Bjellqvist", "Dawson",
		"Lehninger", "Murray"), verbose = FALSE) {
	type = match.arg(type);
	pKa  = aa.pK(type);
	if(is.character(pH)) {
		stop("Did you permute the pH and the peptide sequence?");
	}
	chrg = absoluteCharge.default(x, pH=pH, type = pKa);
	return(chrg);
}



# Note:
# - was initially a Low level function;
#' @export
absoluteCharge = function(x, pH, type) {
	UseMethod("absoluteCharge");
}

# x = String with AA-sequence;
# type = data.frame with pKa values or
#        string name of pKa data-set;
#' @export
absoluteCharge.default = function(x, pH, type) {
	if(is.character(type)) {
		type = match.pK_Names(type);
		pKa  = aa.pK(type);
	} else pKa = type;
	AAs = pKa$AA;
	LEN = length(AAs);
	AAa = AAs[- c(LEN-1, LEN)];
	#
	chrg = sapply(as.seqAA(x), function(x) {
		x = as.tblAA(x, AAs = AAa);
		absoluteCharge.tblAA(x, pH = pH, pKa = pKa);
	});
	return(chrg);
}

# x = Frequency table of charged AAs;
#' @export
absoluteCharge.tblAA = function(x, pH, pKa) {
	charge = aaCharge(pH, pKa);
	charge = sum(charge * x);
	return(charge);
}

aaCharge = function(pH, pKa) {
	div = 1.0 + 10^(pKa$charges * (pH - pKa$pK));
	charge = pKa$charges / div;
}

### Frequency table of AAs:
#' @export
as.tblAA = function(x, AAs = c("C", "D", "E", "H", "K", "R", "Y")) {
	if(inherits(x, "tblAA")) return(x);
	fAA = table(x);
	ids = match(AAs, names(fAA));
	fAA = ifelse(is.na(ids), 0, fAA[ids]);
	fAA = c(fAA, 1, 1); # add: ('c', 'n')
	names(fAA) = c(AAs, "c", "n");
	class(fAA) = c("tblAA", class(x));
	return(fAA);
}

### Convert string to seq of AAs:
#' @export
as.seqAA = function(x) {
	if(inherits(x, "seqAA")) return(x);
	x = strsplit(x, "", fixed = TRUE);
	class(x) = c("seqAA", class(x));
	return(x);
}


### pKa of Amino-Acids
#' @export
aa.pK = function(pKscale = NULL) {
	pKa = list(
	# c = C-terminal; n = N-terminal;
	charges = c(-1, -1, -1, 1, 1, 1, -1, -1, 1),
	EMBOSS = data.frame(
		AA  = c("C", "D", "E", "H", "K", "R", "Y",   "c", "n"),
        pKa = c(8.5, 3.9, 4.1, 6.5, 10.8, 12.5, 10.1, 3.6, 8.6)
	),
	Bjellqvist = data.frame(
		AA  = c("C", "D", "E", "H", "K", "R", "Y",   "c", "n"),
		pKA = c(9, 4.05, 4.45, 5.98, 10, 12, 10,   3.55, 7.50)
	),
	Dawson = data.frame(
		AA  = c("C", "D", "E", "H", "K", "R", "Y",   "c", "n"),
		pKA = c(8.3, 3.9, 4.3, 6.0, 10.5, 12.0, 10.1, 3.2, 8.2)
	),
	Lehninger = data.frame(
		AA  = c("C", "D", "E", "H", "K", "R", "Y",   "c", "n"),
		pKA = c(8.18, 3.65, 4.25, 6, 10.53, 12.48, 10.07,  2.34, 9.69)
	),
	Murray = data.frame(
		AA  = c("C", "D", "E", "H", "K", "R", "Y",   "c", "n"),
		pKA = c(8.33, 3.68, 4.25, 6.0, 11.50, 11.50, 10.07,  2.15, 9.52)
	)
	);
	if(is.null(pKscale)) return(pKa);
	lst = cbind(pKa[[pKscale]], charges = pKa$charges);
	return(lst);
}
match.pK_Names = function(x) {
	nms = c("EMBOSS", "Bjellqvist", "Dawson",
		"Lehninger", "Murray");
	id = pmatch(x, nms);
	nms[id];
}

#################

### Split Mutation into Tokens
# x = "XnnnY", where X, Y = 1 letter codes for AAs;
# - Covers Insertions & Deletions;
#' @export
split.Mutation = function(x, as.upper = TRUE) {
	len = nchar(x);
	if(len < 2) return(NULL);
	tk = regexec("^([A-Za-z])(\\d++)([A-Za-z]*+$)", x, perl = TRUE);
	tk = tk[[1]];
	if(tk[1] == -1) return(NULL);
	nS = tk[-1];
	nE = attr(tk, "match.length")[-1];
	hasAA2 = nE[3] > 0;
	nE = nE + nS - 1;
	nPos = substr(x, nS[2], nE[2]);
	nPos = as.integer(nPos);
	aa1  = substr(x, nS[1], nE[1]);
	aa2  = if(hasAA2) substr(x, nS[3], nE[3]) else "";
	if(as.upper) {
		aa1 = toupper(aa1);
		aa2 = toupper(aa2);
	}
	lst = list(nPos = nPos, Aa1 = aa1, Aa2 = aa2);
	return(lst);
}
# [old]
split.Mutation.old = function(x, as.upper = TRUE) {
	len = nchar(x);
	if(len < 2) return(NULL);
	aa1 = substr(x, 1, 1);
	ch1 = as.integer(charToRaw(aa1));
	if(ch1 < 65 || ch1 > 122) {
		return(NULL); # A, z;
	}
	if(ch1 > 90 && ch1 < 97) {
		return(NULL); # 'Z' < aa1 < 'a';
	}
	if(as.upper) aa1 = toupper(aa1);
	sT2  = substr(x, 2, len);
	sT2  = strsplit(sT2, "(?<=\\d)(?!\\d)", perl=TRUE);
	sT2  = sT2[[1]];
	nPos = sT2[1];
	nPos = as.integer(nPos);
	aa2  = if(length(sT2) == 1) "" else sT2[2];
	if(as.upper) aa2 = toupper(aa2);
	lst  = list(nPos = nPos, Aa1 = aa1, Aa2 = aa2);
	return(lst);
}
