.PHONY: clean

all:\
	data/sorbitolControlDegHg38.rds\
	data/sorbitolControlDegHg19.rds
	
clean:
	rm -rf data/*
	
data/sorbitolControlDegHg38.rds\
data/sorbitolControlDegHg19.rds:\
	rawdata/rna/YAPP_HEK_wt_1_RNApipeSamplesheet.txt\
	rawdata/rna/quant/*/quant.sf\
	scripts/sorbitolControlDeg.R
		mkdir -p data
		Rscript scripts/sorbitolControlDeg.R
