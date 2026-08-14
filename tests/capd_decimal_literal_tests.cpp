#include "../apps/capd_decimal_literal.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace input = hill4bp_capd_input;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    const std::vector<std::string> zero_spellings{
        "0", "-0", "+0.0", "000.000", "0e100", "-0.000e-100"
    };
    for (const std::string& value : zero_spellings) {
        require(input::exact_decimal_zero(value, "test"),
                "zero spelling was not recognized: " + value);
        require(input::exact_decimal_equal(value, "0", "test"),
                "zero spelling was not exactly equal to zero: " + value);
    }

    require(input::exact_decimal_equal("+001.000e0", "1", "test"),
            "exact one normalization failed");
    require(input::exact_decimal_equal("0.003e3", "3", "test"),
            "exact three exponent normalization failed");
    require(input::exact_decimal_equal("30e-1", "3.0", "test"),
            "exact three decimal normalization failed");
    require(!input::exact_decimal_equal(
                "3.0000000000000001", "3", "test"),
            "near-three literal was incorrectly accepted as exact three");
    require(!input::exact_decimal_zero("1e-300", "test"),
            "tiny nonzero literal was incorrectly accepted as exact zero");

    const std::vector<std::string> invalid{
        "", ".", " 0", "0 ", "0x0p0", "nan", "inf", "1e", "--0"
    };
    for (const std::string& value : invalid) {
        bool rejected = false;
        try {
            (void)input::canonical_decimal(value, "test");
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        require(rejected, "invalid literal was accepted: " + value);
    }

    std::cout << "CAPD decimal-literal regression checks passed\n";
    return 0;
}
