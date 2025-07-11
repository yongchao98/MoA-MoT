/*
 * C program for the Wuxing decimal architecture to calculate
 * the percentage of dark matter in the Pandora galaxy.
 */

// #include <wuxing_stdio.h> // Assumed header for frac, printf, etc.
// The provided "frac" type definition is used conceptually.
struct frac {
    signed char n;     // numerator
    unsigned char d;  // denominator
    signed char e;     // exponent
};

// Main function
int main() {
    // === Step 1: Define constants and givens as fractions ===
    // All values are converted to SI units (meters, kg, seconds).
    // The frac type is assumed to be defined in a system header.

    // Gravitational constant G = 6.674e-11 -> approx as 20/3 * 10^-11
    struct frac G = {20, 3, -11};

    // Mass of Sun = 1.989e30 kg -> approx as 2 * 10^30
    struct frac M_sun = {2, 1, 30};

    // Radius r = 10 kpc = 10 * 3.086e19 m -> approx as 3 * 10^20 m
    struct frac r = {3, 1, 20};

    // Velocity v = 200 km/s = 2 * 10^5 m/s
    struct frac v = {2, 1, 5};

    // Luminosity L = 2e9 L_sun
    struct frac lum = {2, 1, 9};

    // Mass-to-light ratio = 3 M_sun/L_sun
    struct frac ml_ratio = {3, 1, 0};

    // === Step 2: Calculate Luminous and Total Mass ===
    // The compiler overloads operators *, /, - for the frac type.
    // The frac library handles simplification to prevent overflow.

    // M_luminous = L * (M/L) * M_sun
    // (2/1e9) * (3/1e0) * (2/1e30) = 12/1e39
    struct frac M_luminous = lum * ml_ratio * M_sun;

    // M_total = (v^2 * r) / G
    // ((2/1e5)^2 * (3/1e20)) / (20/3e-11) = (12/1e30) / (20/3e-11) = 9/5e41
    struct frac M_total = (v * v * r) / G;
    
    // === Step 3: Calculate the dark matter percentage ===

    // frac_dark = 1 - (M_luminous / M_total)
    // 1 - (12/1e39 / 9/5e41) = 1 - (1/15) = 14/15
    struct frac dark_frac = 1 - (M_luminous / M_total);
    
    // percentage = dark_frac * 100
    // (14/15e0) * (100/1e0) is simplified by the library to 14/15e2
    struct frac percentage = dark_frac * 100;
    
    // === Step 4: Convert final fraction to printable decimal value ===
    // We want to print value of frac 'percentage' = (n/d)*10^e rounded to 0.1
    // Calculate (n * 10^(e+1)) / d using integer arithmetic.
    
    int n = percentage.n; // 14
    int d = percentage.d; // 15
    int e = percentage.e; // 2
    
    // Calculate 10^e. Using a loop for general-purpose solution.
    int power_of_10 = 1;
    int i = 0;
    for (i = 0; i < e; i = i + 1) {
        power_of_10 = power_of_10 * 10;
    } // result: 100
    
    // val_x_10 = (14 * 100 * 10) / 15 = 14000 / 15 = 933
    int val_x_10 = (n * power_of_10 * 10) / d; 
    
    int int_part = val_x_10 / 10;   // 93
    int frac_part = val_x_10 % 10; // 3

    // === Step 5: Output the results ===
    // "Output each number in the final equation"
    // The final equation is: P = (1 - M_lum/M_tot) * 100
    printf("Luminous Mass (M_lum): %f kg\n", M_luminous);
    printf("Total Mass (M_tot): %f kg\n", M_total);
    printf("Dark Matter Percentage = (1 - %f / %f) * 100\n", M_luminous, M_total);
    printf("Result: %d.%d%%\n", int_part, frac_part);

    return 0;
}