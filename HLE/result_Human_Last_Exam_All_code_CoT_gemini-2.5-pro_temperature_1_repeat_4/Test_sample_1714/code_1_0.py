import textwrap

# The user wants a C program for the hypothetical Wuxing architecture.
# The request format requires providing this C code within a Python script block.
# This script will print the C code solution.

c_code = """
#include <stdio.h>

/*
 * This program calculates the percentage of dark matter in the Pandora galaxy
 * for the Wuxing decimal-based computer architecture.
 *
 * Wuxing Custom Data Types:
 * - frac: {signed char n, unsigned char d, signed char e} -> (n/d) * 10^e
 * - int: 5 decimal digits (5D)
 * - char: 2 decimal digits (2D)
 *
 * Memory usage for variables is calculated as z.
 * The final answer is given in the format p:z, where p is the percentage.
*/

// Definition for the Wuxing fraction type
struct frac {
    signed char n;
    unsigned char d;
    signed char e;
};

// Main program
int main() {
    /*
     * == Step 1: Define Constants as Fractions ==
     * All constants are approximated to fit into the 'frac' type constraints
     * (n: signed char, d: unsigned char).
    */

    // Gravitational constant G = 6.674e-11, approximated as 20/3 * 10^-11
    struct frac G = {20, 3, -11}; // 6D memory

    // Velocity v = 200 km/s = 2e5 m/s
    struct frac v_ms = {2, 1, 5}; // 6D memory

    // Radius r = 10 kpc = 3.086e20 m, approximated as 3 * 10^20 m
    struct frac r_m = {3, 1, 20}; // 6D memory

    // Luminous Mass = 6e9 * M_sun, with M_sun ~= 2e30 kg. M_lum ~= 12e39 kg
    struct frac M_luminous = {12, 1, 39}; // 6D memory

    /*
     * == Step 2: Calculate the Luminous-to-Total Mass Ratio ==
     * Formula: Percentage = (1 - M_luminous / M_total) * 100
     * where M_total = (v^2 * r) / G.
     * We first calculate the ratio: M_luminous / M_total = (M_luminous * G) / (v^2 * r)
     * This avoids calculating the large total mass directly.
    */

    // Numerator of the ratio: num_ratio = M_luminous * G
    // (12/1e39) * (20/3e-11) = (240/3)e28 = 80/1e28
    struct frac num_ratio = {80, 1, 28}; // 6D memory

    // Denominator of the ratio: den_ratio = v^2 * r
    // v^2 = (2/1e5)^2 = 4/1e10
    struct frac v_sq = {4, 1, 10}; // 6D memory
    // den_ratio = (4/1e10) * (3/1e20) = 12/1e30
    struct frac den_ratio = {12, 1, 30}; // 6D memory

    // The luminous mass ratio. The Wuxing library simplifies this automatically.
    // luminous_ratio = num_ratio / den_ratio = (80/1e28) / (12/1e30) = (80/12)e-2 -> 20/3e-2 -> 1/15e0
    struct frac luminous_ratio = {1, 15, 0}; // 6D memory

    /*
     * == Step 3: Calculate Dark Matter Percentage ==
     * We use integer arithmetic from here to calculate the final percentage value
     * to avoid overflow with the frac type when multiplying by 100.
    */

    // Dark matter ratio = 1 - luminous_ratio = 1 - 1/15 = 14/15
    struct frac dark_ratio = {14, 15, 0}; // 6D memory

    // Extract numerator and denominator for integer calculation
    signed int dr_num = dark_ratio.n;      // 5D memory
    unsigned int dr_den = dark_ratio.d;    // 5D memory

    // Percentage = (dr_num / dr_den) * 100
    signed int p_num = dr_num * 100;       // 14 * 100 = 1400. 5D memory
    signed int p_den = dr_den;             // 15.

    // Calculate integer and first decimal place for printing (e.g., 93.3)
    signed int integer_part = p_num / p_den;      // 1400 / 15 = 93. 5D memory
    signed int remainder = p_num % p_den;         // 1400 % 15 = 5. 5D memory
    unsigned int decimal_part = (remainder * 10) / p_den; // (5 * 10) / 15 = 3. 5D memory

    /*
     * == Step 4: Output the Final Equation and Result ==
     * The final equation is (1 - 1/15) * 100 = 93.3%
    */
    printf("Final Equation: (1 - %d / %d) * 100 = %d.%d%%\\n",
           luminous_ratio.n, luminous_ratio.d, integer_part, decimal_part);

    /*
     * == Step 5: Calculate Total Memory Usage (z) ==
     * 9 frac variables * 6D/frac = 54D
     * 6 int variables * 5D/int = 30D
     * Total z = 54D + 30D = 84D
     */

    return 0;
}
"""
print(textwrap.dedent(c_code).strip())
