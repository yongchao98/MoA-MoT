import textwrap

def generate_wuxing_c_code():
    """
    This function generates the C code for the Wuxing architecture to solve the Pandora galaxy problem.
    """
    c_code = """
    #include <wuxing.h> // Hypothetical header for the Wuxing architecture

    /*
     * This program calculates the percentage of dark matter in the Pandora galaxy
     * on the Wuxing decimal-based architecture.
     *
     * Data:
     * - Velocity (v) = 200 km/s
     * - Radius (r) = 10 kpc
     * - Luminosity = 2e9 L_sun
     * - Mass/Light Ratio = 3 M_sun/L_sun
     *
     * Derived Luminous Mass = (2e9 L_sun) * (3 M_sun/L_sun) = 6e9 M_sun
     *
     * Total Mass Formula: M_total = (v^2 * r) / G
     * The gravitational constant G in astronomical units is approx. 4.30091e-6 kpc*M_sun^-1*(km/s)^2.
     * To keep intermediate calculations within the frac limits (-99 to 99),
     * we approximate G as 13/3 * 10^-6.
     */
    int main() {
        // Variable Definitions & Memory Usage Calculation:
        // 7 frac variables * (2D+2D+2D) = 42D
        // 5 int variables * 5D = 25D
        // Total Memory (z) = 42 + 25 = 67D

        // Use G in units of kpc*M_sun^-1*(km/s)^2, approximated as 13/3 * 10^-6
        frac G = 13/3e-6;

        // Galaxy properties
        frac v = 200;          // Velocity in km/s
        frac r = 10;           // Radius in kpc
        frac M_luminous = 6/1e9; // Luminous mass in solar masses

        // Intermediate variables for calculation
        frac M_total;
        frac lum_ratio;
        frac lum_percent;

        // Equation Step 1: Calculate total mass in solar masses
        // M_total = (v*v*r)/G = (200*200*10)/(13/3e-6) = 400000 / (13/3e-6)
        // This simplifies to 12/13e11 M_sun
        M_total = (v * v * r) / G;

        // Equation Step 2: Calculate luminous matter percentage
        // Luminous % = (M_luminous / M_total) * 100
        // (6e9 / (12/13e11)) * 100 = (6*13/12 * 10^-2) * 100 = 78/12 * 10^-2 * 100 = 13/2
        lum_ratio = M_luminous / M_total;
        lum_percent = lum_ratio * 100; // Result is frac{n=13, d=2, e=0}

        // Deconstruct lum_percent (13/2) to print as a decimal (6.5)
        // and to calculate dark matter percentage without frac overflow.
        int lum_p_int = lum_percent.n / lum_percent.d;     // 13 / 2 = 6
        int lum_p_rem_num = (lum_percent.n % lum_percent.d) * 10; // (1 % 2) * 10 = 10
        int lum_p_frac = lum_p_rem_num / lum_percent.d; // 10 / 2 = 5

        // Equation Step 3: Calculate dark matter percentage from luminous percentage
        // Dark % = 100.0 - Luminous % = 100.0 - 6.5 = 93.5
        // We compute this as (99 - 6) for the integer part and (10 - 5) for the decimal part.
        int dark_p_int = 99 - lum_p_int;  // 99 - 6 = 93
        int dark_p_frac = 10 - lum_p_frac; // 10 - 5 = 5

        // The final equation is 100% - Luminous% = Dark Matter%
        // We output the numbers for the luminous and dark matter percentages.
        printf("Luminous mass percentage is %d.%d percent.\\n", lum_p_int, lum_p_frac);
        printf("The final equation is: 100%% - %d.%d%% = %d.%d%%\\n", lum_p_int, lum_p_frac, dark_p_int, dark_p_frac);
        printf("Percentage of dark matter in Pandora, rounded to 0.1%%, is %d.%d%%.\\n", dark_p_int, dark_p_frac);
        
        return 0;
    }
    """
    print(textwrap.dedent(c_code).strip())

generate_wuxing_c_code()