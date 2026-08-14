#ifndef HILL4BP_CAPD_DECIMAL_LITERAL_HPP
#define HILL4BP_CAPD_DECIMAL_LITERAL_HPP

#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace hill4bp_capd_input {

// Canonical exact value of a finite base-10 literal.  For nonzero values,
// `digits * 10^exponent10` is the unsigned mathematical value; digits has no
// leading or trailing zeroes.  Zero has digits="0", exponent10=0, and no sign.
struct CanonicalDecimal {
    bool negative = false;
    std::string digits = "0";
    long long exponent10 = 0;
};

inline std::runtime_error decimal_error(const std::string& key,
                                        const std::string& text) {
    return std::runtime_error(
        key + " must be one finite base-10 decimal literal (got: " + text + ")");
}

inline CanonicalDecimal canonical_decimal(const std::string& text,
                                          const std::string& key) {
    if (text.empty()) throw decimal_error(key, text);

    std::size_t position = 0;
    bool negative = false;
    if (text[position] == '+' || text[position] == '-') {
        negative = text[position] == '-';
        ++position;
        if (position == text.size()) throw decimal_error(key, text);
    }

    std::string digits;
    bool saw_digit = false;
    while (position < text.size()
           && std::isdigit(static_cast<unsigned char>(text[position]))) {
        saw_digit = true;
        digits.push_back(text[position++]);
    }

    std::size_t fractional_digits = 0;
    if (position < text.size() && text[position] == '.') {
        ++position;
        while (position < text.size()
               && std::isdigit(static_cast<unsigned char>(text[position]))) {
            saw_digit = true;
            digits.push_back(text[position++]);
            ++fractional_digits;
        }
    }
    if (!saw_digit) throw decimal_error(key, text);

    bool exponent_negative = false;
    unsigned long long exponent_magnitude = 0;
    if (position < text.size()
        && (text[position] == 'e' || text[position] == 'E')) {
        ++position;
        if (position < text.size()
            && (text[position] == '+' || text[position] == '-')) {
            exponent_negative = text[position] == '-';
            ++position;
        }
        if (position == text.size()
            || !std::isdigit(static_cast<unsigned char>(text[position]))) {
            throw decimal_error(key, text);
        }
        while (position < text.size()
               && std::isdigit(static_cast<unsigned char>(text[position]))) {
            const unsigned int digit =
                static_cast<unsigned int>(text[position] - '0');
            if (exponent_magnitude
                > (static_cast<unsigned long long>(
                       std::numeric_limits<long long>::max()) - digit) / 10ULL) {
                throw decimal_error(key, text);
            }
            exponent_magnitude = exponent_magnitude * 10ULL + digit;
            ++position;
        }
    }
    if (position != text.size()) throw decimal_error(key, text);

    std::size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw decimal_error(key, text);
    }
    if (consumed != text.size() || !std::isfinite(parsed)) {
        throw decimal_error(key, text);
    }

    const std::size_t first_nonzero = digits.find_first_not_of('0');
    if (first_nonzero == std::string::npos) {
        return CanonicalDecimal{};
    }
    digits.erase(0, first_nonzero);

    if (fractional_digits
        > static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
        throw decimal_error(key, text);
    }
    long long exponent10 = static_cast<long long>(exponent_magnitude);
    if (exponent_negative) exponent10 = -exponent10;
    const long long fractional = static_cast<long long>(fractional_digits);
    if (exponent10 < std::numeric_limits<long long>::min() + fractional) {
        throw decimal_error(key, text);
    }
    exponent10 -= fractional;

    while (digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
        if (exponent10 == std::numeric_limits<long long>::max()) {
            throw decimal_error(key, text);
        }
        ++exponent10;
    }
    return CanonicalDecimal{negative, digits, exponent10};
}

inline bool exact_decimal_equal(const std::string& left,
                                const std::string& right,
                                const std::string& key) {
    const CanonicalDecimal a = canonical_decimal(left, key);
    const CanonicalDecimal b = canonical_decimal(right, key);
    return a.negative == b.negative && a.digits == b.digits
        && a.exponent10 == b.exponent10;
}

inline bool exact_decimal_zero(const std::string& text,
                               const std::string& key) {
    return canonical_decimal(text, key).digits == "0";
}

}  // namespace hill4bp_capd_input

#endif
