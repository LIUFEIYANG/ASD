# automatic process ! DO NOT EDIT !
# setup paths
# ROOT:=$(abspath $(lastword $(MAKEFILE_LIST)))
ROOT:=../..
SRC=$(ROOT)/source
LIB=$(SRC)/lib
APP=$(SRC)/app
BIN=$(ROOT)/bin
OBJ=object
DEP=dependence

# setup targets and libs
lib_include= -I $(LIB)/Eigen -I $(LIB)/gsl/include
#lib_so= -L$(LIB)/gsl/lib -lgsl -lgslcblas -lm
lib_so=
lib_a=
cflags= -Wall
cxxflags= -std=c++11 -Wall -Wno-sign-compare -Wno-unused-variable -Wno-int-in-bool-context
optimazations= -O4 -msse2 -msse3 -msse4
debugoptions= -g

# setup configuration
cc=gcc
cxx=g++

# search targets
targets= $(patsubst %.cpp,%,$(notdir $(wildcard $(APP)/*.cpp))) $(patsubst %.c,%,$(notdir $(wildcard $(APP)/*.c)))

# c/cpp source file name == stored in src/lib & src/app
c_scr_lib=$(wildcard $(LIB)/*.c)
c_src_app=$(wildcard $(APP)/*.c)
cpp_src_lib=$(wildcard $(LIB)/*.cpp)
cpp_src_app=$(wildcard $(APP)/*.cpp)

# c/cpp release/debug object file name == stored in obj/
c_obj_lib_debug=$(patsubst %.c,$(OBJ)/%.debug.o,$(notdir $(c_src_lib)))
c_obj_lib_release=$(patsubst %.c,$(OBJ)/%.o,$(notdir $(c_src_lib)))
c_obj_app_debug=$(patsubst %.c,$(OBJ)/%.debug.o,$(notdir $(c_src_app)))
c_obj_app_release=$(patsubst %.c,$(OBJ)/%.o,$(notdir $(c_src_app)))
cpp_obj_lib_debug=$(patsubst %.cpp,$(OBJ)/%.debug.o,$(notdir $(cpp_src_lib)))
cpp_obj_lib_release=$(patsubst %.cpp,$(OBJ)/%.o,$(notdir $(cpp_src_lib)))
cpp_obj_app_debug=$(patsubst %.cpp,$(OBJ)/%.debug.o,$(notdir $(cpp_src_app)))
cpp_obj_app_release=$(patsubst %.cpp,$(OBJ)/%.o,$(notdir $(cpp_src_app)))

# c/cpp dependence file name == stored in dep/
c_dependence=$(patsubst %.c,$(DEP)/%.d,$(notdir $(c_src_lib)))
c_dependence+=$(patsubst %.c,$(DEP)/%.d,$(notdir $(c_src_app)))
c_dependence+=$(patsubst %.c,$(DEP)/%.debug.d,$(notdir $(c_src_lib)))
c_dependence+=$(patsubst %.c,$(DEP)/%.debug.d,$(notdir $(c_src_app)))
cpp_dependence=$(patsubst %.cpp,$(DEP)/%.d,$(notdir $(cpp_src_lib)))
cpp_dependence+=$(patsubst %.cpp,$(DEP)/%.d,$(notdir $(cpp_src_app)))
cpp_dependence+=$(patsubst %.cpp,$(DEP)/%.debug.d,$(notdir $(cpp_src_lib)))
cpp_dependence+=$(patsubst %.cpp,$(DEP)/%.debug.d,$(notdir $(cpp_src_app)))

# default target == build all the exist application project
.PHONY: all r d clean
all : r d

r : $(patsubst %,%_release,$(notdir $(targets)))

d : $(patsubst %,%_debug,$(notdir $(targets)))

# include dependence file for automatic dependence generation
sinclude $(c_dependence)
sinclude $(cpp_dependence)

# release object
$(OBJ)/%.o:$(LIB)/%.c
	$(cc) $(lib_include) -c $(cflags) $(optimazations) $< -o $@

$(OBJ)/%.o:$(APP)/%.c
	$(cc) $(lib_include) -c $(cflags) $(optimazations) $< -o $@

$(OBJ)/%.o:$(LIB)/%.cpp
	$(cxx) $(lib_include) -c $(cxxflags) $(optimazations) $< -o $@

$(OBJ)/%.o:$(APP)/%.cpp
	$(cxx) $(lib_include) -c $(cxxflags) $(optimazations) $< -o $@

# debug object
$(OBJ)/%.debug.o:$(LIB)/%.c
	$(cc) $(lib_include) -c $(cflags) $(debugoptions) $< -o $@

$(OBJ)/%.debug.o:$(APP)/%.c
	$(cc) $(lib_include) -c $(cflags) $(debugoptions) $< -o $@

$(OBJ)/%.debug.o:$(LIB)/%.cpp
	$(cxx) $(lib_include) -c $(cxxflags) $(debugoptions) $< -o $@

$(OBJ)/%.debug.o:$(APP)/%.cpp
	$(cxx) $(lib_include) -c $(cxxflags) $(debugoptions) $< -o $@

# c dependence file
$(DEP)/%.d:$(LIB)/%.c
	@set -e; rm -f $@ ;\
	$(cc) $(lib_include) -MM $(cflags) $< -MT "$(OBJ)/$(notdir $(<:.c=.o))" -MT "$@" > $@

$(DEP)/%.d:$(APP)/%.c
	@set -e; rm -f $@ ;\
	$(cc) $(lib_include) -MM $(cflags) $< -MT "$(OBJ)/$(notdir $(<:.c=.o))" -MT "$@" > $@

$(DEP)/%.debug.d:$(LIB)/%.c
	@set -e; rm -f $@ ;\
	$(cc) $(lib_include) -MM $(cflags) $< -MT "$(OBJ)/$(notdir $(<:.c=.debug.o))" -MT "$@" > $@

$(DEP)/%.debug.d:$(APP)/%.c
	@set -e; rm -f $@ ;\
	$(cc) $(lib_include) -MM $(cflags) $< -MT "$(OBJ)/$(notdir $(<:.c=.debug.o))" -MT "$@" > $@

# cpp dependence file
$(DEP)/%.d:$(LIB)/%.cpp
	@set -e; rm -f $@ ;\
	$(cxx) $(lib_include) -MM $(cxxflags) $< -MT "$(OBJ)/$(notdir $(<:.cpp=.o))" -MT "$@" > $@

$(DEP)/%.d:$(APP)/%.cpp
	@set -e; rm -f $@ ;\
	$(cxx) $(lib_include) -MM $(cxxflags) $< -MT "$(OBJ)/$(notdir $(<:.cpp=.o))" -MT "$@" > $@

$(DEP)/%.debug.d:$(LIB)/%.cpp
	@set -e; rm -f $@ ;\
	$(cxx) $(lib_include) -MM $(cxxflags) $< -MT "$(OBJ)/$(notdir $(<:.cpp=.debug.o))" -MT "$@" > $@

$(DEP)/%.debug.d:$(APP)/%.cpp
	@set -e; rm -f $@ ;\
	$(cxx) $(lib_include) -MM $(cxxflags) $< -MT "$(OBJ)/$(notdir $(<:.cpp=.debug.o))" -MT "$@" > $@

# targets release
%_release: $(cpp_obj_lib_release) $(OBJ)/%.o
	$(cxx) $^ $(lib_so) $(lib_a) -o $(BIN)/$@

# targets debug
%_debug: $(cpp_obj_lib_debug) $(OBJ)/%.debug.o
	$(cxx) $^ $(lib_so) $(lib_a) -o $(BIN)/$@

# tagets clean
clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*
	rm -f $(DEP)/*
