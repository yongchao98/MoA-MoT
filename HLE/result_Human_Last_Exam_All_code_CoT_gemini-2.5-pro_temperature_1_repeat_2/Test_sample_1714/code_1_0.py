# This Python script generates the C code for the Wuxing platform.
# The user can execute this script to get the required C code.

c_code_for_wuxing = """
#include <stdio.h> // Assuming a standard-like I/O library exists on Wuxing

// The 'frac' type is a special feature of the Wuxing C compiler,
// automatically handling fractional arithmetic. Its definition is:
// struct frac { signed char n; unsigned char d; signed char e; };

int main() {
    // Part 1: Define constants and observed values.
    // All values are stored in the 'frac' type to handle non-integer arithmetic.

    // Velocity v = 200 km/s
    frac v = 200;
    // Radius R = 10 kpc
    frac R = 10;
    // Luminosity L = 2e9 L_sun
    frac Luminosity = 2e9;
    // Mass-to-Light Ratio = 3 (in solar units)
    frac ML_ratio = 3;
    // Gravitational Constant G ≈ 4.3e-6 kpc*(km/s)^2/M_sun.
    // We represent this as 43/10e-6 to fit in the 'frac' type.
    frac G = 43/10e-6;

    // Part 2: Calculate the Luminous and Total Mass of the galaxy.

    // Luminous Mass = Luminosity * Mass-to-Light Ratio
    frac M_luminous = Luminosity * ML_ratio;

    // Total Mass M_total = (v^2 * R) / G
    frac M_total = (v * v * R) / G;

    // Part 3: Calculate the dark matter percentage.
    // Formula: Percentage = (M_total - M_luminous) / M_total * 100
    frac percentage = (M_total - M_luminous) / M_total * 100;

    // Part 4: Print each number in the final equation.
    // The final equation is P = (M_total - M_luminous) / M_total * 100.
    // We print the values for M_total, M_luminous, and the resulting percentage.
    printf("Total Mass (M_sun): %f\\n", M_total);
    printf("Luminous Mass (M_sun): %f\\n", M_luminous);
    printf("Dark Matter Percentage: %f%%\\n", percentage);
    
    return 0;
}
"""

print(c_code_for_wuxing)