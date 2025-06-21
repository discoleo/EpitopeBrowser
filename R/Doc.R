
### Basic Help Page

getHelp = function() {
	HTML("<div>",
		"<p><b>Filter: Rank</b></p>",
		"<p>IEDB recommends a cutoff value of <b>0.5</b> (for HLA-1 epitopes).",
		"A higher cutoff value (up to <b>1 - 1.5</b>) may still work well for epitopes derived from viral structural proteins.",
		"A stringent cutoff value (e.g. <b>0.2</b>) may be only rarely required.</p>",
		"</div>",
		"<div>",
		"<p><b>SubSequence Tab</b></p>",
		"<p>Search all subsequences of a given oligo-peptide for known epitopes.",
		"The epitopes must be loaded as a data set.",
		"The option allEpi.SubSeq specifies which internal data-set is used to search for the epitopes.",
		"By default, epitopes are searched using the unfiltered/full data-set.</p>",
		"</div>"
	);
}