#include <stdio.h>

/*
    This program calculates the percentage of dark matter in the Pandora galaxy
    for the hypothetical Wuxing decimal-based computer architecture.

    Assumptions:
    - The 'frac' type and its associated arithmetic operations (+, -, *, /) are
      pre-defined by the Wuxing C compiler.
    - The compiler can handle conversions from 'frac' to 'int'.
    - The printf function supports the custom '%f' specifier for the 'frac' type.
*/

// The 'frac' struct is defined by the architecture as:
// struct frac {
//     signed char n;     // numerator
//     unsigned char d;  // denominator
//     signed char e;     // exponent
// };
// Value = (n/d) * 10^e

int main() {
    // --- Define Constants as frac types ---
    // The values are approximated to fit within 'char' limits for n and d.

    // Velocity v = 200 km/s = 2e5 m/s
    // frac{n=2, d=1, e=5}
    frac v = 2/1e5;

    // Radius r = 10 kpc = 3.086e20 m. Approximated as 3.1e20 m.
    // To fit n in signed char, this is represented as 31 * 10^19.
    // frac{n=31, d=1, e=19}
    frac r = 31/1e19;

    // Gravitational Constant G = 6.674e-11. Approximated as 6.7e-11.
    // Represented as (67/10) * 10^-11.
    // frac{n=67, d=10, e=-11}
    frac G = 67/10e-11;

    // Luminosity L = 2e9 L_sun
    // frac{n=2, d=1, e=9}
    frac L = 2/1e9;

    // Mass-to-Light Ratio = 3
    // frac{n=3, d=1, e=0}
    frac M_L_ratio = 3/1e0;

    // Solar Mass M_sun = 1.989e30 kg. Approximated as 2e30 kg.
    // frac{n=2, d=1, e=30}
    frac M_sun = 2/1e30;

    // --- Calculation Steps ---

    // 1. Calculate Total Mass: M_total = (v^2 * r) / G
    frac v_squared = v * v;
    frac M_total = (v_squared * r) / G;

    // 2. Calculate Luminous Mass: M_luminous = L * M_L_ratio * M_sun
    frac M_luminous = L * M_L_ratio * M_sun;

    // 3. Calculate Dark Matter Mass: M_dark = M_total - M_luminous
    frac M_dark = M_total - M_luminous;

    // 4. Calculate Dark Matter Percentage: (M_dark / M_total) * 100
    frac hundred = 100/1e0;
    frac percentage_raw = (M_dark / M_total) * hundred;

    // 5. Round the result to one decimal place using integer/frac arithmetic.
    // Method: floor(value * 10 + 0.5)
    frac ten = 10/1e0;
    frac half = 5/10e0; // Represents 0.5

    frac temp_for_rounding = (percentage_raw * ten) + half;
    int rounded_val_times_10 = (int)temp_for_rounding; // Cast to int truncates (floors)

    int final_integer_part = rounded_val_times_10 / 10;
    int final_decimal_part = rounded_val_times_10 % 10;

    // --- Final Output ---
    // As requested, output the numbers in the final equation.
    printf("Dark Matter Percentage Equation:\n");
    printf("(Dark Mass / Total Mass) * 100 = Percentage\n");
    printf("(%f / %f) * 100 = %d.%d%%\n", M_dark, M_total, final_integer_part, final_decimal_part);

    return 0;
}
<<<93.3:105>>>