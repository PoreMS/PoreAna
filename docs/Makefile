# Minimal makefile for Sphinx documentation
#

SPHINXOPTS    =
ifdef OS
   SPHINXBUILD   = python -msphinx
else
   SPHINXBUILD   = python3 -msphinx
endif
SPHINXPROJ    = PoreAna
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
