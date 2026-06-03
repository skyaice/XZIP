# Compatibility wrapper for the historical Debug build entry.
# The canonical build configuration lives in the repository-root Makefile.

.PHONY: all clean

all:
	$(MAKE) -C .. all BUILD_DIR=Debug/build

clean:
	$(MAKE) -C .. clean BUILD_DIR=Debug/build
