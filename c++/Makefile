CC=g++
GTEST_DIR=/Users/usuyama/codes/class_test/hapmuc/tests/libs/gtest-1.7.0
CFLAGS=-c -Wall -I./ -I$(GTEST_DIR)/include
LDFLAGS=-L$(GTEST_DIR) -lgtest
SOURCES := test.cpp Read.cpp VariantFromSAM.cpp
OBJECTS=$(SOURCES:.cpp=.o)
test: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
