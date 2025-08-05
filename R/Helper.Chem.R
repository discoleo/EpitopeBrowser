
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

