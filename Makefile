# © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
# Switzerland, Laboratory of Experimental Biophysics, 2017
# See the LICENSE.txt file for more details.

CC := g++
SRCDIR := PolymerCpp/core
BUILDDIR := build
TARGETDIR := lib
TARGET := $(TARGETDIR)/PolymerCpp.so

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -std=c++11 -O2 -fPIC
LIB := $(shell python-config --ldflags) -fopenmp
INC := -I include/ $(shell python-config --cflags)

$(TARGET): $(OBJECTS)
	@mkdir -p $(TARGETDIR)
	@echo " Linking..."
	@echo " $(CC) -shared $^ -o $(TARGET) $(LIB)"; $(CC) -shared $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo $(LIB)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $< $(LIB)"; $(CC) $(CFLAGS) $(INC) -c -o $@ $< $(LIB)

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
