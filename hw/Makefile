# Alternative to CMake. Usage: make [CXX=g++]; make test
CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
BUILD := build

all: $(BUILD)/uf_decode $(BUILD)/test_uf

$(BUILD)/uf_decode: src/main.cpp src/uf.cpp src/uf.hpp
	mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp src/uf.cpp

$(BUILD)/test_uf: tests/test_uf.cpp src/uf.cpp src/uf.hpp
	mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ tests/test_uf.cpp src/uf.cpp

test: $(BUILD)/test_uf
	./$(BUILD)/test_uf

clean:
	rm -rf $(BUILD)

.PHONY: all test clean
