

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
