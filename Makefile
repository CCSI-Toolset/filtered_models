# A dispatch makefile for the filtered models project running make in
# each of it's component subdirs
SUBDIRS := gas-particle-cylinder_2d_models \
           gas-solid-3d-models
CLEAN_TARGETS := $(SUBDIRS:%=clean-%)

# Where Jenkins should checkout ^/projects/common/trunk/
COMMON := .ccsi_common
LEGAL_DOCS := LEGAL \
           CCSI_TE_LICENSE.txt

.PHONY: all $(SUBDIRS) $(CLEAN_TARGETS)

all: $(LEGAL_DOCS) $(SUBDIRS)

$(LEGAL_DOCS):
	@if [ -d $(COMMON) ]; then \
	  cp $(COMMON)/$@ .; \
	else \
	  svn -q export ^/projects/common/trunk/$@; \
	fi

$(SUBDIRS):
	@$(MAKE) -sC $@

clean: $(CLEAN_TARGETS)
	@rm -rf $(LEGAL_DOCS)

$(CLEAN_TARGETS):
	@$(MAKE) -sC $(@:clean-%=%) clean
