
# Source:
# 1. http://www.allelefrequencies.net/top10freqs.asp
# 2. Sacchi N, Castagnetta M, Miotti V, Garbarino L, Gallina A.
#    High-resolution analysis of the HLA-A, -B, -C and -DRB1 alleles and national and
#    regional haplotype frequencies based on 120,926 volunteers from the Italian
#    Bone Marrow Donor Registry.
#    HLA. 2019 Sep;94(3):285-295. https://doi.org/10.1111/tan.13613. PMID: 31207125;
# 3. Eberhard HP, Schmidt AH, Mytilineos J, Fleischhauer K, Müller CR.
#    Common and well-documented HLA alleles of German stem cell donors by
#    haplotype frequency estimation.
#    HLA. 2018 Oct;92(4):206-214. https://doi.org/10.1111/tan.13378. PMID: 30117303;
#
#' @export
hla = function() {
	list(
	# for Europe:
	A = data.frame(
		HLA = c("A*01:01", "A*02:01", "A*03:01", "A*11:01", "A*24:02", "A*26:01", "A*32:01",
			# Italy: A68, A30, A02 {14,26,03};
			# Germany: A02 {05, 06}, A*33;
			"A*68:01", "A*30:01", "A*02:14", "A*02:26", "A*02:03",
			"A*02:05", "A*02:06", "A*33:01", "A*33:03"),
		Type = "A",
		Freq = c(0.124, 0.2624, 0.1185, 0.059, 0.101, 0.0385, 0.0385,
			0.029, 0.027, 0.000095, 0.000083, 0.000062,
			0.0037, 0.00052, 0.0006, 0.00019)),
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
