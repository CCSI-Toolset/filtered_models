# A simple makefile for creating the Sorbentfit compressed tarball file
VERSION := 2016.02.00
PRODUCT := FilteredModels_GasParticleCylinder2D
LICENSE := CCSI_TE_LICENSE_$(PRODUCT).txt
PKG_DIR := CCSI_$(PRODUCT)_$(VERSION)
PACKAGE := $(PKG_DIR).tgz

PARENT_PAYLOAD := ../LEGAL
PARENT_PAYLOAD_LOCAL := $(shell echo $(PARENT_PAYLOAD) | sed 's:../::g')

PAYLOAD := $(wildcard docs/*.pdf) \
	   model/     \
	   example1   \
	   example2   \
	   test/      \
	   $(LICENSE)

# Get just the top part (not dirname) of each entry so cp -r does the right thing
PAYLOAD_TOPS := $(foreach v,$(PAYLOAD),$(shell echo $v | cut -d'/' -f1))
# And the payload with the PKG_DIR prepended
PKG_PAYLOAD := $(addprefix $(PKG_DIR)/, $(PAYLOAD) $(PARENT_PAYLOAD_LOCAL))

# OS detection & changes
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
  MD5BIN := md5sum
endif
ifeq ($(UNAME), Darwin)
  MD5BIN := md5
endif
ifeq ($(UNAME), FreeBSD)
  MD5BIN := md5
endif

.PHONY: all clean

all: $(PACKAGE)

# Make compressed tar file without timestamp (gzip -n) so md5sum
# doesn't change if the payload hasn't
$(PACKAGE): $(PAYLOAD)
	@mkdir $(PKG_DIR)
	@cp $(PARENT_PAYLOAD) $(PKG_DIR)
	@cp -r $(PAYLOAD_TOPS) $(PKG_DIR)
	@tar -cf - $(PKG_PAYLOAD) | gzip -n > $(PACKAGE)
	@$(MD5BIN) $(PACKAGE)
	@rm -rf $(PKG_DIR) $(LICENSE)

$(LICENSE): ../CCSI_TE_LICENSE.txt 
	@sed "s/\[SOFTWARE NAME \& VERSION\]/$(PRODUCT) v.$(VERSION)/" < ../CCSI_TE_LICENSE.txt > $(LICENSE)

clean:
	@rm -rf $(PACKAGE) $(PKG_DIR) $(LICENSE)
