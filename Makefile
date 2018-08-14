CXX=g++

DIR_IMPL=Implementation
DIR_TEST=Test
DIR_EXAMPLE=Example
DIR_LIB_TEST=Test/Lib
DIR_BIN=Bin
DIR_BIN_EXAMPLE=$(DIR_BIN)/Example
DIR_BIN_TEST=$(DIR_BIN)/Test
INCLUDES=-I $(DIR_IMPL)
INCLUDES_TEST=-I $(DIR_IMPL) -I $(DIR_LIB_TEST)
CPPFLAGS=-Wall -Wextra -std=c++11 -O2 -flto $(INCLUDES)
CPPFLAGS_TEST=-Wall -Wextra -std=c++11 -g -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer $(INCLUDES_TEST)
CPPMACROS=
CPPMACROS_TEST=-DCNRC_DEBUG_AVL_SEQUENCE
CPPMACROS_ENABLE_HDT=-DCNRC_ENABLE_HDT

HPPFILES_IMPL=$(wildcard $(DIR_IMPL)/*.hpp)

BINARIES_EXAMPLE=$(DIR_BIN_EXAMPLE)/DecrementalGC \
                 $(DIR_BIN_EXAMPLE)/DecrementalGC-HDT \
                 $(DIR_BIN_EXAMPLE)/DecrementalGMCC \
                 $(DIR_BIN_EXAMPLE)/DecrementalGMCC-HDT

BINARIES_TEST=$(DIR_BIN_TEST)/AVLSequenceTest \
              $(DIR_BIN_TEST)/EulerTourTreeTest \
              $(DIR_BIN_TEST)/EulerTourTreeSpanningForestTest \
              $(DIR_BIN_TEST)/HDTSpanningForestTest \
			  $(DIR_BIN_TEST)/DecrementalMCCTest \
			  $(DIR_BIN_TEST)/DecrementalMCCTest-HDT

usage :
	@echo "Please specify a target : 'make example' or 'make test'."

example : $(BINARIES_EXAMPLE)

$(DIR_BIN_EXAMPLE)/DecrementalGMCC : $(DIR_EXAMPLE)/DecrementalGMCC.cpp $(HPPFILES_IMPL) | $(DIR_BIN_EXAMPLE)
	$(CXX) $(CPPFLAGS) $(CPPMACROS) -o $@ $<

$(DIR_BIN_EXAMPLE)/DecrementalGMCC-HDT : $(DIR_EXAMPLE)/DecrementalGMCC.cpp $(HPPFILES_IMPL) | $(DIR_BIN_EXAMPLE)
	$(CXX) $(CPPFLAGS) $(CPPMACROS) $(CPPMACROS_ENABLE_HDT) -o $@ $<

$(DIR_BIN_EXAMPLE)/DecrementalGC : $(DIR_EXAMPLE)/DecrementalGC.cpp $(HPPFILES_IMPL) | $(DIR_BIN_EXAMPLE)
	$(CXX) $(CPPFLAGS) $(CPPMACROS) -o $@ $<

$(DIR_BIN_EXAMPLE)/DecrementalGC-HDT : $(DIR_EXAMPLE)/DecrementalGC.cpp $(HPPFILES_IMPL) | $(DIR_BIN_EXAMPLE)
	$(CXX) $(CPPFLAGS) $(CPPMACROS) $(CPPMACROS_ENABLE_HDT) -o $@ $<

test : $(BINARIES_TEST)

$(DIR_BIN_TEST)/AVLSequenceTest : $(DIR_TEST)/AVLSequenceTest.cpp $(DIR_IMPL)/AVLSequence.hpp | $(DIR_BIN_TEST)
	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) -o $@ $<
	@echo "** AVLSequenceTest may take some time.**"
	$@

$(DIR_BIN_TEST)/EulerTourTreeTest : $(DIR_TEST)/EulerTourTreeTest.cpp $(DIR_BIN_TEST)/AVLSequenceTest $(DIR_IMPL)/EulerTourTree.hpp | $(DIR_BIN_TEST)
	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) -o $@ $<
	$@

$(DIR_BIN_TEST)/EulerTourTreeSpanningForestTest : $(DIR_TEST)/EulerTourTreeSpanningForestTest.cpp $(DIR_BIN_TEST)/EulerTourTreeTest $(DIR_IMPL)/EulerTourTreeSpanningForest.hpp | $(DIR_BIN_TEST)
	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) -o $@ $<
	$@

$(DIR_BIN_TEST)/HDTSpanningForestTest : $(DIR_TEST)/HDTSpanningForestTest.cpp $(DIR_BIN_TEST)/EulerTourTreeTest $(DIR_IMPL)/HDTSpanningForest.hpp | $(DIR_BIN_TEST)
	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) -o $@ $<
	$@

$(DIR_BIN_TEST)/DecrementalMCCTest : $(DIR_TEST)/DecrementalMCCTest.cpp $(DIR_BIN_TEST)/EulerTourTreeSpanningForestTest $(DIR_IMPL)/DecrementalMCCModule.hpp $(DIR_IMPL)/DecrementalMCC.hpp | $(DIR_BIN_TEST)
	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) -o $@ $<
	$@

$(DIR_BIN_TEST)/DecrementalMCCTest-HDT : $(DIR_TEST)/DecrementalMCCTest.cpp $(DIR_BIN_TEST)/HDTSpanningForestTest $(DIR_IMPL)/DecrementalMCCModule.hpp $(DIR_IMPL)/DecrementalMCC.hpp | $(DIR_BIN_TEST)
	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) $(CPPMACROS_ENABLE_HDT) -o $@ $<
	$@

#$(BINARIES_TEST) : % : %.cpp $(HPPFILES_IMPL)
#	$(CXX) $(CPPFLAGS_TEST) $(CPPMACROS_TEST) -o $@ $@.cpp
#	$@

Bin/Example :
	mkdir -p Bin/Example

Bin/Test :
	mkdir -p Bin/Test

clean :
	rm -f $(BINARIES_EXAMPLE) $(BINARIES_TEST)

.PHONY : example test clean usage
