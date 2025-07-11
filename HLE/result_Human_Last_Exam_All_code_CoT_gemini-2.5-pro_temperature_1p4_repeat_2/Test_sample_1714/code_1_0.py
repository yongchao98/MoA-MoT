import textwrap

# This Python script generates the C code for the Wuxing architecture.
# The C code calculates the percentage of dark matter in the Pandora galaxy.

c_code = """
#include <stdio.h>
// Fictional header for the Wuxing architecture, providing the 'frac' type.
#include <wuxing.h>

int main() {
    // This program calculates the dark matter percentage of Pandora galaxy.
    // All values and calculations use the special 'frac' data type.
    
    // ---- Constants and Initial Data ----
    // G = Gravitational Constant in (km/s)^2 * kpc / M_sun. ~4.3e-6 is written as 43/10 * 10^-7.
    frac G = 43/10e-7;
    // r = Radius in kpc.
    frac r = 10;
    // v = Velocity in km/s.
    frac v = 200;
    // L_galaxy = Luminosity in L_sun units.
    frac L_galaxy = 2e9;
    // M_L_ratio = Mass-to-light ratio in M_sun/L_sun units.
    frac M_L_ratio = 3;

    // ---- Calculation ----
    // To ensure numerical stability and optimize memory, the calculation is done
    // using the formula: Percentage = (1 - LuminousFraction) * 100.

    // 1. Calculate Luminous Mass Fraction = (Luminous Mass) / (Total Mass)
    // M_luminous = L_galaxy * M_L_ratio
    // M_total = (v^2 * r) / G
    frac luminous_fraction = (L_galaxy * M_L_ratio) / ((v * v * r) / G);

    // 2. Calculate the final percentage.
    // The integer literals 1 and 100 are automatically converted to 'frac' type.
    frac percentage = (1 - luminous_fraction) * 100;

    // ---- Output ----
    // To meet the requirement, we output the numbers that form the final equation.
    // The final step is: Percentage = (1 - luminous_fraction) * 100
    printf("The final calculation is: (1 - Luminous_Fraction) * 100\\n");
    printf("(1 - %f) * 100 = %f %%\\n", luminous_fraction, percentage);

    return 0;
}
"""

print(textwrap.dedent(c_code).strip())