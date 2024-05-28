library(ggplot2)

default_pointsize <- 8
quartzFonts(sans = quartzFont(c("Helvetica Neue Light", "Helvetica Neue", "Helvetica Light","Helvetica")))

ggp2_magic_number <- ggplot2::.pt * 72.27/96
gtext_magic_number <- 0.36

if( capabilities()['aqua'] ) {
  default_font <- "Helvetica Neue Light"
} else {
  default_font <- "Helvetica"
}

pdf_out <- function( file,
	dpi=NA_real_, family=default_font, pointsize = default_pointsize, ... ){
	if( capabilities()['aqua'] ){
		# preferred device type with configurable fonts
		quartz( type="pdf", file=file, dpi=dpi, family=family, pointsize,  ...)
		par( family=family )

	} else {
		# fallback device type with default Helvetica font
		pdf( file, family=family, pointsize, ...)
	}
}
