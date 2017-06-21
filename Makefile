# A dispatch makefile for the filtered models project running make in
# each of it's component subdirs
SUBDIRS := gas-particle-cylinder_2d_models \
           gas-solid-3d-models
CLEAN_TARGETS := $(SUBDIRS:%=clean-%)

.PHONY: all $(SUBDIRS) $(CLEAN_TARGETS)

all: $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) -sC $@

clean: $(CLEAN_TARGETS)

$(CLEAN_TARGETS):
	@$(MAKE) -sC $(@:clean-%=%) clean
