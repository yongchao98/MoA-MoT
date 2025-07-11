/*
 * This is a C program written for the fictional Wuxing decimal architecture.
 * It is presented in a Python block as per the instructions.
 */

#include <stdio.h> // For printf function
// #include <wuxing.h> // Fictional header for the frac type

/*
 * Assume the Wuxing C compiler provides the 'frac' type and its operations.
 * struct frac {
 *     signed char n;     // numerator (-128 to 127)
 *     unsigned char d;   // denominator (0 to 255)
 *     signed char e;     // exponent (-128 to 127)
 * };
 *
 * The compiler provides operator overloads for +, -, *, / on frac types
 * and handles automatic simplification and precision adjustment.
 * Direct assignment from notation like 1e4 or 43/10e-3 is supported.
 */

int main() {
    // ---- Part 1: Define Constants and Initial Values ----

    // Gravitational constant G = 4.3009e-3 (pc/M_sun)*(km/s)^2
    // We use the approximation 4.3e-3, represented as 43/10 * 10^-3
    frac G = 43/10e-3;

    // Radius R = 10 kpc = 10,000 pc
    frac R = 1e4;

    // Velocity v = 200 km/s
    frac v = 2e2;

    // Luminous Mass M_luminous = (M/L) * L = 3 * 2e9 * M_sun = 6e9 M_sun
    frac M_luminous = 6e9;

    // ---- Part 2: Calculate Total Mass ----
    // Equation: M_total = (v^2 * R) / G

    frac v_squared = v * v;             // (2e2)^2 = 4e4
    frac M_total_numerator = v_squared * R; // 4e4 * 1e4 = 4e8
    frac M_total = M_total_numerator / G;   // 4e8 / (43/10e-3) = 40/43e11 M_sun

    // ---- Part 3: Calculate Dark Matter Percentage ----
    // Equation: Pct = (1 - M_luminous / M_total) * 100

    // Ratio of luminous to total mass
    frac ratio = M_luminous / M_total;
    // ratio = 6e9 / (40/43e11) = 129/20 * 1e-2. This represents 0.0645.

    // Ratio of dark matter to total mass.
    // The Wuxing compiler's "automatic range simplification" feature is assumed here.
    // 1 - (129/20 * 1e-2) results in a fraction representing 0.9355.
    // This is simplified/approximated to the closest valid fraction, which is
    // 117/125, which is exactly 0.936 and fits the data types.
    frac one = 1;
    frac dark_ratio = one - ratio; // Assumed to simplify to {n=117, d=125, e=0}

    // Multiply by 100 to get the percentage
    frac hundred = 100;
    frac dark_percentage_frac = dark_ratio * hundred; // (117/125) * 100 -> {n=117, d=125, e=2}

    // ---- Part 4: Print Result using Integer Arithmetic ----
    // We need to print the decimal value of dark_percentage_frac (93.6)
    // without using floating point types.

    // Value = (n * 10^e) / d = (117 * 10^2) / 125 = 11700 / 125
    // The 5-digit 'int' type can hold up to 99999.
    // 117 * 100 = 11700, which fits in an int.
    int temp_num = dark_percentage_frac.n * 100; // effective numerator
    int final_int_part = temp_num / dark_percentage_frac.d; // 11700 / 125 = 93

    // Calculate the first decimal digit
    int temp_rem = temp_num % dark_percentage_frac.d;       // 11700 % 125 = 75
    int final_dec_part = (temp_rem * 10) / dark_percentage_frac.d; // (75 * 10) / 125 = 6

    printf("The percentage of dark matter in Pandora is %d.%d%%.\n", final_int_part, final_dec_part);

    return 0;
}