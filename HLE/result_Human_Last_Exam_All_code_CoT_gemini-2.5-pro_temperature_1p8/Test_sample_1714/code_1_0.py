import textwrap

def generate_c_code():
    """
    Generates the C code for the Wuxing architecture to solve the dark matter problem.
    """
    c_code = """\
    #include <stdio.h>

    /*
     * This program calculates the percentage of dark matter in the Pandora galaxy
     * using the Wuxing decimal-based architecture.
     *
     * We assume the Wuxing C compiler provides the following:
     * - A 'frac' type: struct frac { signed char n; unsigned char d; signed char e; };
     * - Automatic conversion of integers to 'frac' (e.g., 100 becomes 1/1e2).
     * - Overloaded arithmetic operators (+, -, *, /) for 'frac' types.
     * - A custom printf specifier '%f' that displays a frac, e.g., "93/1e9".
     */

    // Forward declaration of our custom printing function
    void print_frac_as_decimal(struct frac f);

    int main() {
        // --- GIVENS AND CONSTANTS ---
        // Mass calculation constant K = 2.325e5 solar_mass*s^2/(kpc*km^2)
        // Represented as 93/40 * 10^5 to fit in 'frac' type constraints.
        frac K = 93/40e5;

        // Radius r = 10 kpc
        frac r = 10;

        // Velocity v = 200 km/s, represented as 2 * 10^2
        frac v = 2/1e2;
        
        // Luminous Mass M_lum = 6e9 solar masses
        frac M_lum = 6/1e9;

        // --- CALCULATION VARIABLES ---
        frac M_total;
        frac v_sq;
        frac dark_mass;
        frac dark_ratio;
        frac hundred;
        frac dark_percent;

        // --- COMPUTATION ---
        // 1. Calculate Total Mass: M_total = K * r * v^2
        // v_sq will be 4/1e4. M_total will simplify to 93/1e9 (9.3e10 M_sun)
        v_sq = v * v;
        M_total = K * r * v_sq;

        // 2. Calculate Dark Mass: dark_mass = M_total - M_lum
        // Result is 87/1e9
        dark_mass = M_total - M_lum;
        
        // 3. Calculate dark matter percentage
        // dark_ratio will simplify from (87/1e9)/(93/1e9) to 29/31
        dark_ratio = dark_mass / M_total;
        
        // The Wuxing compiler converts the integer 100 to a frac
        hundred = 100;
        dark_percent = dark_ratio * hundred; // Result is 29/31e2

        // --- OUTPUT ---
        // The prompt requires printing the numbers in the final equation.
        // We use the simplified numerators from the frac objects for clarity.
        // M_total = 93/1e9, dark_mass = 87/1e9.
        printf("Percentage = (%d / %d) * 100 = ", dark_mass.n, M_total.n);
        
        // Print the final percentage value, rounded to 0.1%
        print_frac_as_decimal(dark_percent);
        printf("\\n");
        
        return 0;
    }

    /**
     * @brief Prints a frac value as a decimal number rounded to 0.1 precision.
     * This function uses integer arithmetic exclusively to avoid overflow on the
     * 5-digit 'int' type by performing division in a loop.
     * @param f The fraction to print.
     */
    void print_frac_as_decimal(struct frac f) {
        // --- LOCAL VARIABLES ---
        int int_part = 0;
        int rem;
        int dec1, dec2;
        char sign = 1;
        int i;

        // Extract components to avoid modifying the original 'frac'
        signed char num = f.n;
        unsigned char den = f.d;
        signed char exp = f.e;

        if (num < 0) {
            sign = -1;
            num = -num;
        }

        // Calculate the integer part before the decimal point.
        // Example: For 29/31e2, this loop calculates 93.
        int_part = num / den;
        rem = num % den;
        for (i = 0; i < exp; i++) {
            int_part = int_part * 10;
            rem = rem * 10;
            int_part += rem / den;
            rem = rem % den;
        }

        // Calculate first decimal digit
        rem *= 10;
        dec1 = rem / den;
        rem %= den;

        // Calculate second decimal digit for rounding
        rem *= 10;
        dec2 = rem / den;

        // Perform rounding
        if (dec2 >= 5) {
            dec1++;
        }
        
        // Handle carry-over if rounding pushes the decimal part to 10
        if (dec1 >= 10) {
            int_part++;
            dec1 = 0;
        }
        
        if (sign == -1) {
            printf("-");
        }
        printf("%d.%d%%", int_part, dec1);
    }
    """
    # Print the C code block as requested
    print("```c\n" + textwrap.dedent(c_code) + "```")

generate_c_code()