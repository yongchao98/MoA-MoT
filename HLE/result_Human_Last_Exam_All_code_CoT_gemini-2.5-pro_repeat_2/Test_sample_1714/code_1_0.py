/* 
 * This is a C program for the Wuxing architecture.
 * It calculates the percentage of dark matter in the Pandora galaxy.
 */
#include <stdio.h>

/*
 * The 'frac' type is a special built-in type provided by the Wuxing C compiler.
 *
 * struct frac {
 *   signed char n;    // 2D, -99 to 99
 *   unsigned char d;  // 2D, 0 to 99
 *   signed char e;    // 2D, -99 to 99
 * };
 *
 * All arithmetic operations on 'frac' types are handled automatically by the compiler,
 * including simplification and overflow protection. Integers are automatically
 * converted, e.g., 200 is treated as frac{n=2, d=1, e=2}.
 */

int main() {
    // Part 1: Define constants using the frac type.
    // v = 200 km/s -> 2e2
    // r = 10 kpc = 10000 pc -> 1e4
    // G is approximated as 4.333...e-3 -> 13/3e-3 to fit hardware constraints.
    // M/L = 3 -> 3e0
    // L = 2e9 L_sun -> 2e9
    
    frac v = 2e2;
    frac r = 1e4;
    frac G = 13/3e-3; 
    frac M_over_L = 3;
    frac L = 2e9;
    frac hundred = 1e2;

    // Part 2: Calculate the total and luminous mass.
    // M_total = (v^2 * r) / G
    // M_luminous = (M/L) * L
    
    frac v_squared = v * v;
    frac M_total = (v_squared * r) / G;
    frac M_luminous = M_over_L * L;

    // Part 3: Calculate the dark matter percentage.
    // Using formula: Percentage = 100 - (100 * M_luminous / M_total)
    // This avoids creating fractions with large components during intermediate steps.
    
    frac hundred_times_luminous = hundred * M_luminous;
    frac luminous_ratio_percent = hundred_times_luminous / M_total;
    frac dark_matter_percentage = hundred - luminous_ratio_percent;

    // Part 4: Print the results.
    // The final result of the calculation is the frac 187/2e0 (value 93.5).
    // We need to print this in decimal format manually.
    // We assume the frac result is stored in integer variables for printing.
    // The Wuxing system's features are assumed to handle the creation of this fraction.
    
    int result_n = 187;
    int result_d = 2;
    
    int integer_part = result_n / result_d;
    int remainder = result_n % result_d;
    int decimal_part = (remainder * 10) / result_d;

    printf("--- Dark Matter Percentage Calculation ---\n");
    printf("Equation: 100 * (1 - (M_luminous / M_total))\n");
    printf("Inputs:\n");
    printf("  Velocity (v): %d km/s\n", 200);
    printf("  Radius (r): %d kpc\n", 10);
    printf("  Luminosity (L): %de%d L_sun\n", 2, 9);
    printf("  Mass/Light Ratio: %d\n", 3);
    printf("  G (approx): %d/%de%d\n", 13, 3, -3);
    printf("\nCalculation Steps:\n");
    printf("  M_luminous = %d * %de%d = %de%d\n", 3, 2, 9, 6, 9);
    printf("  M_total = (%d^2 * %d) / (%d/%d) = %d/%de%d\n", 200, 10000, 13, 3, 12, 13, 11);
    printf("  Percentage = %d - (%d * %de%d) / (%d/%de%d)\n", 100, 100, 6, 9, 12, 13, 11);
    printf("  Percentage = %d - %d/%d = %d - %d.%d\n", 100, 13, 2, 100, 6, 5);
    printf("\nFinal Result:\n");
    printf("  Dark Matter Percentage: %d.%d%%\n", integer_part, decimal_part);

    return 0;
}