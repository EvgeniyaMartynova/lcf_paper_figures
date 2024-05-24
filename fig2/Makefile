.SECONDARY:
.DELETE_ON_ERROR:

all : figure2.pdf

figure2.pdf : figure2.tex images_witness.txt 
	latexmk -use-make -xelatex $<


images_witness.txt : figure2.R ../settings.R ../utils_lcf.R
	Rscript $< && touch images_witness.txt