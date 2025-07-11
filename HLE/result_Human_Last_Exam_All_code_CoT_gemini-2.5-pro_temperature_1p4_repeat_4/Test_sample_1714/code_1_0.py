// This is a C program for the Wuxing architecture.
#include <stdio.h>

// Definition of the special fraction type for Wuxing
// This would typically be in a system header file.
struct frac {
    signed char n;     // numerator
    unsigned char d;  // denominator
    signed char e;     // exponent
};

// Assume the compiler provides built-in operators (+, -, *, /) for the frac type
// that handle automatic simplification and overflow protection.

int main() {
    // --- Variable Declarations and Memory Calculation ---
    // Total memory usage (z) will be calculated by summing the size of all variables.
    // 1 frac = 3 * char = 3 * 2D = 6D.
    // 1 int = 5D.
    // We declare 9 frac variables and 6 int variables.
    // z = (9 * 6D) + (6 * 5D) = 54D + 30D = 84D.

    // Part 1: Define constants as frac types.
    // Selections are made to keep numerators/denominators small to avoid overflow.
    struct frac v = {2, 1, 5};                  // Velocity: 200 km/s = 2e5 m/s
    struct frac r = {31, 10, 20};               // Radius: 10 kpc ~= 3.1e20 m
    struct frac G = {2, 3, -10};                // Gravitational Constant: G ~= 6.67e-11 ~= 2/3 e-10
    struct frac lum_factor = {6, 1, 9};         // Luminous mass factor from L and M/L ratio
    struct frac M_sun = {2, 1, 30};             // Mass of the Sun: M_sun ~= 2e30 kg
    
    // Part 2: Print the equation and values used.
    printf("Calculating Dark Matter Percentage = (1 - Luminous_Mass / Total_Mass) * 100\n");
    printf("Using the following frac approximations:\n");
    printf("Velocity (v)      : %d/%d * 10^%d m/s\n", v.n, v.d, v.e);
    printf("Radius (r)        : %d/%d * 10^%d m\n", r.n, r.d, r.e);
    printf("G Constant        : %d/%d * 10^%d\n", G.n, G.d, G.e);
    printf("Luminous Factor   : %d/%d\n", lum_factor.n * M_sun.n, lum_factor.d * M_sun.d);
    printf("Luminous Exponent : 10^%d * 10^%d\n", lum_factor.e, M_sun.e);
    printf("\n");

    // Part 3: Calculate the luminous mass to total mass ratio.
    // The order of operations is critical to prevent intermediate overflow.
    struct frac temp_calc;
    temp_calc = lum_factor * M_sun; // (6/1e9) * (2/1e30) = 12/1e39
    temp_calc = temp_calc / v;      // (12/1e39) / (2/1e5) = 6/1e34
    temp_calc = temp_calc / v;      // (6/1e34)  / (2/1e5) = 3/1e29
    temp_calc = temp_calc / r;      // (3/1e29)  / (31/10e20) = 30/31e9
    struct frac lum_ratio = temp_calc * G; // (30/31e9) * (2/3e-10) = 20/31e-1 = 2/31

    // Part 4: Calculate the final percentage.
    struct frac one = {1, 1, 0};
    struct frac hundred = {1, 1, 2};
    struct frac dark_ratio = one - lum_ratio; // 1 - 2/31 = 29/31
    struct frac dark_percentage_frac = dark_ratio * hundred; // (29/31) * 100

    // Part 5: Convert frac to a rounded decimal and print.
    // We want to print (29/31)*100 rounded to one decimal place.
    // Use integer arithmetic to perform rounding.
    int num = dark_percentage_frac.n; // 29
    int den = dark_percentage_frac.d; // 31
    // The value is num/den * 10^e, which is 29/31 * 100
    
    // To round to 1 decimal place, calculate value*100, add 5, then divide by 10.
    // (29 * 100 * 100 + 31/2) / 31 / 10
    int result_times_100 = num * 100 * 100 / den; // (290000 / 31) = 9354
    int rounded_result_times_10 = (result_times_100 + 5) / 10; // (9354 + 5) / 10 = 935
    
    int integer_part = rounded_result_times_10 / 10; // 93
    int decimal_part = rounded_result_times_10 % 10; // 5

    printf("--- Calculation Result ---\n");
    printf("Dark matter percentage: %d.%d%%\n", integer_part, decimal_part);

    return 0;
}
<<<93.5:84>>>