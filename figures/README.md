# Figures

This folder contains all the scripts and files required to generate the figures. The data figures were generated using `R` and saved as separate panels, before being manually assembled into the figures presented using either `Inkscape` or `LibreOffice Impress`.

To generate all figure panels, run:
`R -e rmarkdown::render('figures.Rmd')`

To compile main figures PNG images into a PDF, run:
`convert figures/*.png figures.pdf`

To convert supplement.odp to a PDF, run:
`soffice --headless --convert-to pdf figures/supplement.odp`
