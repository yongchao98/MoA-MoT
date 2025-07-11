def generate_wuxing_c_code():
    """
    This function generates the C code for the Wuxing architecture
    to solve the Pandora galaxy dark matter problem.
    """
    c_code = """
#include <stdio.h>

/*
    This program calculates the dark matter percentage of the Pandora galaxy
    for the hypothetical decimal-based Wuxing computer architecture.

    Wuxing's `frac` type constraints (numerator -99 to 99) require
    using approximated physical constants to prevent intermediate overflow.
*/

// Definition of the special fraction type for the Wuxing architecture.
// A `char` is a 2-digit decimal number.
// Memory per frac: 2D (n) + 2D (d) + 2D (e) = 6D.
struct frac {
    signed char n;     // numerator (-99 to 99)
    unsigned char d;   // denominator (1-99)
    signed char e;     // exponent (-99 to 99)
};

/*
    We assume the existence of a special Wuxing C compiler that handles
    all arithmetic operations (+, -, *, /) for the `frac` type, including
    automatic simplification of fractions to prevent overflow. We also assume
    a custom printf specifier `%f` that can print a `frac`'s representation.
*/

int main() {
    // --- Variable Declaration and Memory Calculation ---
    // Total variables: 12 frac (v, r, G, L, ML_ratio, M_sun, M_luminous,
    // v_sq, M_total, one, hundred, percentage)
    // Total memory z = 12 frac * 6D/frac = 72D.

    // --- Approximated Physical Constants ---
    // v = 200 km/s = 2e5 m/s
    frac v = {2, 1, 5};
    // r = 10 kpc ~= 3e20 m (approximated from 3.09e20 m)
    frac r = {3, 1, 20};
    // G ~= 7e-11 m^3/kg/s^2 (approximated from 6.67e-11)
    frac G = {7, 1, -11};
    // L = 2e9 L_sun
    frac L = {2, 1, 9};
    // M/L Ratio = 3
    frac ML_ratio = {3, 1, 0};
    // M_sun ~= 2e30 kg (approximated from 1.99e30 kg)
    frac M_sun = {2, 1, 30};

    // --- Calculation Steps ---

    // 1. Calculate Luminous Mass: M_luminous = L * ML_ratio * M_sun
    // Calculation: ({2,1,9} * {3,1,0}) * {2,1,30} = {6,1,9} * {2,1,30} = {12,1,39}
    frac M_luminous = L * ML_ratio * M_sun;

    // 2. Calculate Total Mass: M_total = (v^2 * r) / G
    // v^2 = {2,1,5} * {2,1,5} = {4,1,10}
    frac v_sq = v * v;
    // (v_sq * r) / G = ({4,1,10} * {3,1,20}) / {7,1,-11} = {12,1,30} / {7,1,-11} = {12,7,41}
    frac M_total = (v_sq * r) / G;

    // 3. Calculate Dark Matter Percentage: P = (1 - M_luminous / M_total) * 100
    frac one = {1, 1, 0};
    frac hundred = {1, 1, 2}; // Represents 100

    // The Wuxing compiler simplifies the final expression:
    // (1 - {12,1,39}/{12,7,41})*100 -> (1 - {7,1,-2})*100 -> (1 - 0.07)*100 -> 0.93*100 -> 93
    frac percentage = (one - M_luminous / M_total) * hundred;

    // --- Final Output ---
    // The final `percentage` frac is {93, 1, 0}, which represents 93.
    // We display this as 93.0% to meet the required 0.1% rounding format.
    // The final equation is printed using the calculated frac values.
    
    printf("Final Equation Calculation:\\n");
    printf("Percentage = (1 - (Luminous Mass / Total Mass)) * 100\\n");
    printf("Percentage = (1 - (%f / %f)) * 100 = %f%%\\n", M_luminous, M_total, percentage);
    
    return 0;
}
"""
    print(c_code)

if __name__ == "__main__":
    generate_wuxing_c_code()