.SECONDARY:
.DELETE_ON_ERROR:

all : figure1.pdf

figure1.pdf : figure1.tex images_witness.txt 
	latexmk -use-make -xelatex $<


images_witness.txt : figure1.R ../settings.R ../utils_lcf.R
	Rscript $< && touch images_witness.txt