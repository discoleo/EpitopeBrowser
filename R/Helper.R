

### Basic Functions

### Strings

### Text => Values of Range
splitRange = function(x) {
	txt = strsplit(x, "[-, /]+")[[1]];
	val = round(as.numeric(txt));
	val = val[! is.na(val)];
	return(val);
}

# Used for MHC-2: Core
#' @export
findPos = function(what, x, fixed = TRUE) {
	npos = sapply(seq_along(what), function(id) {
		npos = regexpr(what[id], x[id], fixed = fixed);
		attributes(npos) = NULL;
		return(npos);
	})
	return(npos);
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

# MHC-2 vs MHC-1:
isMHC2 = function(x) {
	idMHC2 = which(grepl("netmhciipan_el.score", names(x)));
	isMHC2 = length(idMHC2) > 0;
	return(isMHC2);
}
getCore = function(x) {
	x$netmhciipan_el.core;
}

### Read Epitopes
# sep = separator used in the csv file;
read.epi = function(file, hla.strip = TRUE, sep = ",") {
	x = read.csv(file, header = TRUE, sep = sep);
	if(is.numeric(x[,1])) {
		# IEDB: new format / unedited;
		hasMultiSeq = length(unique(x[,1])) > 1;
		# MHC-2 vs MHC-1:
		idMHC2 = which(grepl("netmhciipan_el.score", names(x)));
		isMHC2 = length(idMHC2) > 0;
		#
		if(isMHC2) corePP = getCore(x);
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
		# MHC-2 Core:
		if(isMHC2) x$Cs = findPos(corePP, x$Peptide);
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

### HLA Regions


#' @export
merge.HLARegions3 = function(
		x = list(Hu = "Hungary"),
		y = list(De = "De", It = "Italy"), ...) {
	x   = c(y, x);
	hla = merge.HLARegions(x, y = NULL, ...);
	return(hla);
}
#' @export
merge.HLARegions = function(
		x = list(De = "De", It = "Italy", Hu = "Hungary"),
		y = NULL, ...) {
	regions = x;
	nms = c("HLA", "Freq");
	sfx = paste0("Freq.", names(regions));
	x = hla(regions[[1]]);
	y = hla(regions[[2]])[, nms];
	names(x)[3] = sfx[1];
	names(y)[2] = sfx[2];
	x = merge(x, y, by = "HLA", all = TRUE);
	# Region 3:
	y = hla(regions[[3]])[, nms];
	y = y[! is.na(y$HLA), ];
	# Exclude Level 3:
	isL3 = grepl("\\:[0-9]++\\:[0-9]++$", y$HLA, perl = TRUE);
	isL2 = grepl("\\:[0-9]++$", y$HLA, perl = TRUE);
	isL1 = TRUE;
	if(any(isL3) || all(isL2)) {
		isL1 = FALSE;
		x$HLA.L1 = x$HLA; # hack: use L2;
		y$HLA[isL3] = sub("\\:[0-9]++$", "", y$HLA[isL3], perl = TRUE);
	} else {
		# Level 1:
		x$HLA.L1 = sub("\\:[0-9]++$", "", x$HLA, perl = TRUE);
		# Exclude Level 2:
		# - inferred L2 frequencies;
		y = y[! grepl("\\:[0-9]++$", y$HLA, perl = TRUE), ];
	}
	names(y) = c("HLA.L1", sfx[3]);
	x = merge(x, y, by = "HLA.L1", all = TRUE);
	# if(! isL1) x$HLA.L1 = sub("\\:[0-9]++$", "", x$HLA.L1, perl = TRUE);
	print(str(x))
	return(x);
}

#####################

### Epitopes: Summary

names.quant = function(x) {
	names(x)[grepl("^Rq", names(x))];
}

summary.epi = function(x, isMultiSeq = NULL, probs = c(0, 0.5, 1)) {
	if(is.null(isMultiSeq)) {
		isMultiSeq = match("Seq", names(x));
		isMultiSeq = ! is.na(isMultiSeq);
	}
	# Unique Set:
	colNms = if(isMultiSeq) "Seq" else NULL;
	colNms = c("Peptide", colNms);
	datUnq   = unique(x[colNms]);
	# Stats:
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
