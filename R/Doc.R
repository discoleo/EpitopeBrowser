
### Basic Help Page

getHelp = function() {
	HTML(paste0("<div>",
		"<p><b>Rank</b></p>",
		"<p>IEDB recommends a cutoff value of <b>0.5</b> (for HLA-1 epitopes).
		A higher cutoff value (up to <b>1 - 1.5</b>) may still work well for epitopes derived from viral structural proteins. A stringent cutoff value (e.g. <b>0.2</b>) may be only rarely required.</p>",
		"</div>"
	));
}