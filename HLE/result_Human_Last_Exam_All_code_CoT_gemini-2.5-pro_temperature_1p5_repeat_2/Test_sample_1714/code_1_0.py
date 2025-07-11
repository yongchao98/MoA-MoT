def generate_wuxing_c_code():
    """
    This function generates the C code for the Wuxing computer architecture
    to calculate the dark matter percentage of the Pandora galaxy.
    """
    c_code = """#include <stdio.h>

// This program calculates the percentage of dark matter in the Pandora galaxy
// on the hypothetical Wuxing decimal-based computer architecture.
// All calculations are performed using the special 'frac' data type.

int main() {
    // Step 1: Define constants based on the problem description.
    // The Wuxing 'frac' type is defined as {signed char n, unsigned char d, signed char e}.
    // To avoid numerator overflow (max 127), the gravitational constant G is
    // approximated from 4.3e-6 to 4e-6.

    // v = 200 km/s  -> 2/1 * 10^2
    /%frac v = 2/1e2;

    // r = 10 kpc -> 1/1 * 10^1
    /%frac r = 1/1e1;

    // G ≈ 4e-6 (km/s)^2 * kpc / M_solar -> 4/1 * 10^-6
    /%frac G = 4/1e-6;

    // Luminous Mass (M_lum) = (2e9 L_sun) * (3 M_sun/L_sun) = 6e9 M_sun -> 6/1 * 10^9
    /%frac m_lum = 6/1e9;
    
    // Constant for percentage calculation -> 100/1 * 10^0 or 1/1 * 10^2
    /%frac hundred = 1/1e2;


    // Step 2: Calculate the Total Mass: M_total = (v^2 * r) / G
    // The Wuxing compiler handles the fraction arithmetic automatically.
    // Result: M_total = (4/1e5) / (4/1e-6) = 1/1e11
    /%frac m_total = (v * v * r) / G;


    // Step 3: Calculate the Luminous Matter Percentage: (M_lum / M_total) * 100
    // Result: ((6/1e9) / (1/1e11)) * 100 = (6/1e-2) * 100 = 6/1e0, which is 6.
    /%frac lum_perc = (m_lum / m_total) * hundred;
    

    // Step 4: Calculate the Dark Matter Percentage: 100 - Luminous %
    // Result: 100 - 6 = 94.
    /%frac dark_perc = hundred - lum_perc;


    // Step 5: Print the final equation and result rounded to one decimal place.
    // Cast frac results to Wuxing 'int' type for printing.
    int h_int = (int)hundred;
    int l_int = (int)lum_perc;
    int d_int = (int)dark_perc;

    printf("Dark Matter Calculation:\\n");
    printf("%d%% (Total) - %d%% (Luminous) = %d.0%% (Dark Matter)\\n", h_int, l_int, d_int);

    return 0;
}
"""
    print(c_code)

if __name__ == "__main__":
    generate_wuxing_c_code()