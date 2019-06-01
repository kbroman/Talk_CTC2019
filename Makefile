STEM = ctc2019
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

FIGS = Figs/hs.pdf \
	   Figs/overview.pdf \
	   Figs/do_genome.pdf \
	   Figs/geno_reconstruct.pdf \
	   Figs/mapped_reads.pdf \
	   Figs/dist_matrix.pdf \
	   Figs/min_v_self.pdf \
	   Figs/detailed_mixups.pdf \
	   Figs/lrt_v_contam.pdf \
	   Figs/lrt_v_contam_detail.pdf \
	   Figs/lrt_v_contam_expand.pdf

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $<

Figs/%.pdf: R/%.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

notes: $(STEM)_withnotes.pdf
all: $(STEM).pdf notes web

$(STEM)_withnotes.tex: $(STEM).tex Ruby/createVersionWithNotes.rb
	Ruby/createVersionWithNotes.rb $(STEM).tex $(STEM)_withnotes.tex

$(STEM)_withnotes.pdf: $(STEM)_withnotes.tex header.tex $(FIGS)
	xelatex $(STEM)_withnotes
	pdfnup $(STEM)_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv $(STEM)_withnotes-nup.pdf $(STEM)_withnotes.pdf

web: $(STEM).pdf notes
	scp $(STEM)_withnotes.pdf adhara.biostat.wisc.edu:Website/presentations/$(STEM)_withnotes.pdf
	scp $(STEM).pdf adhara.biostat.wisc.edu:Website/presentations/$(STEM).pdf
