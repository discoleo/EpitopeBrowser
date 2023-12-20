
# Source:
# http://www.allelefrequencies.net/top10freqs.asp
#
#' @export
hla = function() {
	list(
	# for Europe:
	A = data.frame(
		HLA = c("A*01:01", "A*02:01", "A*03:01", "A*11:01", "A*24:02", "A*26:01", "A*32:01"),
		Type = "A",
		Freq = c(0.124, 0.2624, 0.1185, 0.059, 0.101, 0.0385, 0.0385)),
	B = data.frame(
		HLA = c("B*07:02", "B*08:01", "B*15:01", "B*18:01", "B*35:01", "B*44:02", "B*44:03", "B*51:01"),
		Type = "B",
		Freq = c(0.0796, 0.075, 0.0465, 0.0549, 0.06, 0.06, 0.0418, 0.075)),
	C = data.frame(
		HLA = c("C*02:09", "C*04:01", "C*05:01", "C*06:02", "C*07:01", "C*07:02", "C*12:03"),
		Type = "C",
		Freq = c(0.0693, 0.129, 0.0702, 0.0906, 0.144, 0.111, 0.0652)
	)
	)
}

#' @export
as.data.frame.hla = function() {
	hla = hla();
	rbind(hla$A, hla$B, hla$C);
}
