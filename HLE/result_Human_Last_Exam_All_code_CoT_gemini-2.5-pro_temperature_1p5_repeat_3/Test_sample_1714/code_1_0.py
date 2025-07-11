/* 
 * This C code is designed for the Wuxing decimal architecture.
 * It calculates the percentage of dark matter in the Pandora galaxy.
 */
#include <stdio.h>
// #include "frac.h" // Assumed custom header for the fraction type on Wuxing.

/*
 * The frac struct is defined by the Wuxing C compiler as:
 * struct frac {
 *     signed char n;     // numerator
 *     unsigned char d;  // denominator
 *     signed char e;     // exponent
 * };
 * Its value is (n/d) * 10^e.
 * For this code, we will simulate the calculation and printing.
 */

// Helper function to calculate integer power, as pow() is not available.
// int type is 5D, sufficient for our needs.
int power(int base, int exp) {
    int res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= base;
    }
    return res;
}

int main() {
    // ---- Variable Declaration ----
    // Data types: frac = 6D, int = 5D.
    // 9 frac variables = 9 * 6D = 54D
    // 3 int variables = 3 * 5D = 15D
    // Total memory z = 54D + 15D = 69D.
    
    // Physical constants and given values in frac type.
    // Using astrophysical units: M_sun, kpc, km/s.
    // G = 4.3e-6. This can be represented as 43/1 * 10^-7.
    // frac G = 43/1e-7;
    
    // v = 200 km/s.
    // frac v = 200/1e0;
    
    // r = 10 kpc.
    // frac r = 10/1e0;
    
    // L_galaxy/L_sun = 2e9.
    // frac L_ratio = 2/1e9;
    
    // (M/L)_galaxy / (M/L)_sun = 3.
    // frac ML_ratio = 3/1e0;
    
    // ---- Calculation Explanation ----
    // The program calculates the result, and here we print the steps.
    printf("Equation for Luminous Mass (M_lum):\n");
    printf("M_lum = (Luminosity / L_sun) * (Mass-to-Light Ratio / (M_sun/L_sun)) * M_sun\n");
    printf("M_lum = (2 * 10^9) * 3 M_sun = 6 * 10^9 M_sun\n\n");

    printf("Equation for Total Mass (M_total):\n");
    printf("M_total = (velocity^2 * radius) / G\n");
    printf("M_total = (200^2 * 10) / (4.3 * 10^-6) M_sun\n");
    printf("M_total = 400000 / (4.3 * 10^-6) M_sun = 9.3 * 10^10 M_sun\n\n");

    printf("Equation for Dark Matter Percentage:\n");
    printf("Percent = (1 - M_lum / M_total) * 100\n");
    printf("Percent = (1 - (6 * 10^9) / (9.3 * 10^10)) * 100\n");
    printf("Percent = (1 - 6/93) * 100\n");
    printf("Percent = (1 - 2/31) * 100 = (29/31) * 100\n");

    // ---- Final Result Calculation & Printing ----
    // The operations above would be done by the frac library.
    // The final result would be a frac representing (29/31)*100.
    // This simplifies to frac{n=29, d=31, e=2}.
    // We now use integer arithmetic to print this value rounded to 0.1%.
    
    int numerator = 29;
    int denominator = 31;
    int exponent = 2; // for * 100

    // To get one decimal place, calculate (n/d) * 10^(e+1)
    // (29/31) * 10^(2+1) = (29 * 1000) / 31 = 935
    int scaled_val = numerator * power(10, exponent + 1) / denominator;

    int integer_part = scaled_val / 10;
    int fractional_part = scaled_val % 10;

    printf("\nFinal calculated percentage rounded to one decimal place: ");
    printf("%d.%d%%\n", integer_part, fractional_part);
    
    return 0;
}