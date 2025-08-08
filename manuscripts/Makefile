MD_FILES := overview_ITC_LLMs.md cNMA_methodology.md MLNMR_methodology.md NMI_methodology.md MAIC_methodology.md STC_methodology.md
PDF_FILES := $(MD_FILES:.md=.pdf)

PANDOC := pandoc
CSL := ieee.csl
TEMPLATE := eisvogel

all: pdf

pdf: $(PDF_FILES)

%.pdf: %.md
	$(PANDOC) --from markdown+tex_math_dollars+escaped_line_breaks \
	  --pdf-engine=xelatex \
	  --filter pandoc-crossref \
	  -V geometry:margin=1in \
	  -V fontsize=11pt \
	  -V linestretch=1.15 \
	  --citeproc --csl=$(CSL) \
	  -o $@ $<

clean:
	rm -f $(PDF_FILES)
