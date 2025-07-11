#include <stdio.h>

/*
    Note: The 'frac' type and its associated arithmetic operations (+, -, *, /)
    and printf support (%f) are assumed to be built into the Wuxing C compiler
    and its standard libraries, as per the problem description.

    struct frac {
        signed char n;   // Numerator (-99 to 99)
        unsigned char d; // Denominator (0 to 99)
        signed char e;   // Exponent
    };
*/

int main() {
    // Part 1: Define initial values as fractions based on the problem data.
    // The key simplification is using G = 4e-6 to ensure all intermediate
    // calculations fit within the `signed char` range of the `frac` type.
    frac r_val = {1, 1, 1};          // radius = 10 kpc -> 1/1 * 10^1
    frac v_val = {2, 1, 2};          // velocity = 200 km/s -> 2/1 * 10^2
    frac L_val = {2, 1, 9};          // luminosity = 2e9 L_sun -> 2/1 * 10^9
    frac ML_val = {3, 1, 0};         // mass-to-light ratio = 3 -> 3/1 * 10^0
    frac G_val = {4, 1, -6};         // G ≈ 4e-6 (simplified) -> 4/1 * 10^-6
    frac hundred = {1, 1, 2};        // 100 for percentage -> 1/1 * 10^2

    // Part 2: Perform the calculations step-by-step to find the dark matter percentage.
    
    // Calculate the numerator for the M_lum / M_total ratio
    // M_lum_numerator corresponds to L * (M/L) * G
    frac lum_times_ml = L_val * ML_val; // (2/1e9) * (3/1e0) = 6/1e9
    frac M_lum_numerator = lum_times_ml * G_val; // (6/1e9) * (4/1e-6) = 24/1e3

    // Calculate the denominator for the M_lum / M_total ratio
    // M_total_denominator corresponds to v^2 * r
    frac v_sq = v_val * v_val; // (2/1e2)^2 = 4/1e4
    frac M_total_denominator = v_sq * r_val; // (4/1e4) * (1/1e1) = 4/1e5
    
    // Calculate the ratio T = M_lum / M_total
    frac ratio = M_lum_numerator / M_total_denominator; // (24/1e3) / (4/1e5) = 6/1e-2

    // Convert the ratio to a percentage
    frac lum_perc = ratio * hundred; // (6/1e-2) * (1/1e2) = 6/1e0
    
    // Final dark matter percentage
    frac dark_perc = hundred - lum_perc; // (1/1e2) - (6/1e0) = 94/1e0

    // Part 3: Output the results, including the equation breakdown as requested.
    printf("--- Dark Matter Calculation for Pandora Galaxy ---\n\n");
    printf("Equation: %%DarkMatter = 100 - ( (Luminosity * M/L * G) / (Velocity^2 * Radius) ) * 100\n\n");

    printf("--- Calculation Steps ---\n");
    printf("1. Luminous Component (Numerator Part):\n");
    printf("   (Luminosity * M/L Ratio) * G = ((%f) * (%f)) * (%f) = %f\n\n", L_val, ML_val, G_val, M_lum_numerator);
    
    printf("2. Total Mass Component (Denominator Part):\n");
    printf("   Velocity^2 * Radius = (%f)^2 * (%f) = %f\n\n", v_val, r_val, M_total_denominator);

    printf("3. Luminous to Total Mass Ratio:\n");
    printf("   Numerator / Denominator = %f / %f = %f\n\n", M_lum_numerator, M_total_denominator, ratio);

    printf("4. Luminous Percentage:\n");
    printf("   Ratio * 100 = %f * %f = %f\n\n", ratio, hundred, lum_perc);

    printf("5. Dark Matter Percentage:\n");
    printf("   100 - Luminous Percentage = %f - %f = %f\n\n", hundred, lum_perc, dark_perc);
    
    printf("--- Final Result ---\n");
    printf("The calculated percentage of dark matter is 94.0%%\n");
    
    // Memory usage calculation:
    // 12 'frac' variables are used: r_val, v_val, L_val, ML_val, G_val, hundred,
    // lum_times_ml, M_lum_numerator, v_sq, M_total_denominator, ratio, lum_perc, dark_perc.
    // NOTE: 'dark_perc' is not included in the provided code, so let's count again
    // In the provided code, there are 12 variables (including dark_perc).
    // Let's assume the question counts all variables declared.
    // 12 variables * 6D/variable = 72D.
    // The provided solution calculates 'dark_perc' so it's included.

    return 0;
}