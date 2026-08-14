CXX ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude
LDFLAGS ?=

BUILD_DIR := build
LIB_OBJECTS := \
	$(BUILD_DIR)/ode.o \
	$(BUILD_DIR)/crtbp.o \
	$(BUILD_DIR)/hill4bp.o

.PHONY: all capd clean test

all: $(BUILD_DIR)/hill4bp_poincare

# Requires CAPD_CONFIG=/path/to/CAPD/build/bin/capd-config (or capd-config on
# PATH).  This builds both the dense reconnaissance and interval-proof apps.
capd:
	BUILD_DIR="$(abspath $(BUILD_DIR))" bash apps/build_capd_apps.sh

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/ode.o: src/ode.cpp include/crtbp_cpp/ode.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/crtbp.o: src/crtbp.cpp include/crtbp_cpp/crtbp.hpp include/crtbp_cpp/ode.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/hill4bp.o: src/hill4bp.cpp include/crtbp_cpp/hill4bp.hpp include/crtbp_cpp/ode.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/hill4bp_poincare.o: apps/hill4bp_poincare.cpp include/crtbp_cpp/hill4bp.hpp include/crtbp_cpp/ode.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/hill4bp_poincare: $(LIB_OBJECTS) $(BUILD_DIR)/hill4bp_poincare.o
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/hill4bp_regression_tests.o: tests/hill4bp_regression_tests.cpp include/crtbp_cpp/hill4bp.hpp include/crtbp_cpp/ode.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/hill4bp_regression_tests: $(LIB_OBJECTS) $(BUILD_DIR)/hill4bp_regression_tests.o
	$(CXX) $(LDFLAGS) $^ -o $@

test: $(BUILD_DIR)/hill4bp_regression_tests $(BUILD_DIR)/capd_decimal_literal_tests
	$(BUILD_DIR)/hill4bp_regression_tests
	$(BUILD_DIR)/capd_decimal_literal_tests
	bash tests/test_neck_sweep.sh
	bash tests/test_professional_neck_study.sh

$(BUILD_DIR)/capd_decimal_literal_tests: tests/capd_decimal_literal_tests.cpp apps/capd_decimal_literal.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(BUILD_DIR)
