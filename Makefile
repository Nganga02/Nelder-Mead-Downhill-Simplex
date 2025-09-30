CC = gcc
CFLAGS = -Wall -std=c17 -pedantic -O2 -I$(INCDIR)
LDFLAGS = -lm

#dir
INCDIR = include
SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES = demo-2D.c\
	xmalloc.c \
	nelder-mead.c

SOURCES02 = demo-energy.c\
	xmalloc.c\
	nelder-mead.c

SRCFILES = 	$(addprefix $(SRCDIR)/,$(SOURCES))

DEMOSRCFILES = $(addprefix $(SRCDIR)/,$(SOURCES02))

OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCFILES))

DEMOOBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(DEMOSRCFILES))

HEADERS = $(wildcard $(INCDIR)/*.h)

TARGET01 =  $(BINDIR)/demo-2D

TARGET02 = $(BINDIR)/demo-energy


demo-2D: $(TARGET01)

demo-energy: $(TARGET02)
# Make directories if not exist
$(OBJDIR):
	@mkdir -p $@

$(BINDIR):
	@mkdir -p $@

# Link the executable
$(TARGET01) : $(OBJECTS) | $(BINDIR)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# For the demo truss files
$(TARGET02) : $(DEMOOBJS) | $(BINDIR)
	$(CC) $(DEMOOBJS) -o $@ $(LDFLAGS)

# Compile source files into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -g -c $< -o $@

# Clean build
clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean