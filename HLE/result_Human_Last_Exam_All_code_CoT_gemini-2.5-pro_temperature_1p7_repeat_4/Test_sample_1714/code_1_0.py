import math

def generate_wuxing_c_code():
    """
    This function generates the C code for the fictional Wuxing computer
    and calculates the final percentage and memory usage.
    """

    # --- Part 1: Calculation based on Wuxing constraints ---

    # Approximations that fit within the 'frac' type's 'signed char' numerator (-128 to 127)
    # v = 2e5, L = 2e9, M/L = 3 are used directly.
    # G ≈ 20/3 * 10^-11 -> frac G = 20/3e-11
    # R ≈ 3e20 -> frac R = 3e20
    # M_sun ≈ 2e30 -> frac M_sun = 2e30

    # These lead to a ratio M_lum/M_total of 1/15.
    # The percentage is (1 - 1/15) * 100 = (14/15) * 100 = 93.333...%
    
    p = round((14 / 15) * 100, 1)

    # --- Part 2: Calculate Memory Usage (z) ---
    
    # 13 frac variables: G, v, R, M_sun, L, M_L_ratio, one, hundred,
    # M_total, M_luminous, lum_total_ratio, dm_frac, dm_perc
    num_frac_vars = 13
    size_of_frac = 6 # 2D for n + 2D for d + 2D for e
    mem_frac = num_frac_vars * size_of_frac

    # 7 int variables for printing logic:
    # i, multiplier, int_part, rem1, dec_part1, rem2, dec_part2
    num_int_vars = 7
    size_of_int = 5 # 5D per int
    mem_int = num_int_vars * size_of_int

    z = mem_frac + mem_int
    
    # --- Part 3: Generate the C code as a string ---

    c_code = """
#include <wuxing.h> // Assumed standard library for Wuxing

// This program calculates the percentage of dark matter in the Pandora galaxy
// using the constraints of the Wuxing decimal architecture.

int main() {
    // Part 1: Define constants and initial values using 'frac' type
    // Approximations are chosen to fit within 'frac' numerator/denominator limits.
    frac v = 2e5;             // Velocity: 200 km/s = 2e5 m/s
    frac R = 3e20;            // Radius: 10 kpc approx 3e20 m
    frac G = 20/3e-11;        // Gravitational Constant: approx 6.67e-11
    frac L = 2e9;             // Luminosity: 2e9 L_sun
    frac M_L_ratio = 3;       // Mass-to-Light ratio
    frac M_sun = 2e30;        // Solar Mass: approx 2e30 kg
    
    // Part 2: Perform the calculations
    // The Wuxing C compiler's frac library handles the complex arithmetic automatically.
    
    // M_total = v^2 * R / G
    frac M_total = (v * v * R) / G; // Result: 9/5e41

    // M_luminous = L * M/L_ratio * M_sun
    frac M_luminous = L * M_L_ratio * M_sun; // Result: 12/1e39

    // Percentage = (1 - M_luminous / M_total) * 100
    frac one = 1;
    frac lum_total_ratio = M_luminous / M_total; // Result: 1/15e0
    frac dm_frac = one - lum_total_ratio;        // Result: 14/15e0
    
    frac hundred = 1e2;
    frac dm_perc = dm_frac * hundred; // Result: 14/15e2

    // Part 3: Print the results step-by-step and the final rounded percentage
    printf("--- Dark Matter Calculation Steps ---\\n");
    printf("1. Luminous Mass (M_lum) = %f kg\\n", M_luminous);
    printf("2. Total Mass (M_total) = %f kg\\n", M_total);
    printf("3. Luminous/Total Ratio = %f\\n", lum_total_ratio);
    printf("4. Dark Matter Fraction = %f\\n", dm_frac);
    printf("5. Dark Matter Percentage (unformatted) = %f %%\\n", dm_perc);
    printf("-------------------------------------\\n");
    
    // Part 4: Manually format and print the final result rounded to 0.1%
    // Note: Integer types are 5D, so they can hold values up to 49999.
    // Calculations must not overflow this limit. dm_perc = 14/15e2
    
    int int_part = (dm_perc.n * 100) / dm_perc.d;     // (14*100)/15 = 93
    int rem1 = (dm_perc.n * 100) % dm_perc.d;       // 1400 % 15 = 5
    int dec_part1 = (rem1 * 10) / dm_perc.d;          // (5*10)/15 = 3
    int rem2 = (rem1 * 10) % dm_perc.d;             // 50 % 15 = 5
    int dec_part2 = (rem2 * 10) / dm_perc.d;          // (5*10)/15 = 3

    // Round up the first decimal place if the second is 5 or greater
    if (dec_part2 >= 5) {
        dec_part1++;
        if (dec_part1 == 10) { // Handle cascading round, e.g., 93.99 -> 94.0
            dec_part1 = 0;
            int_part++;
        }
    }

    printf("Final Rounded Result:\\n");
    printf("The percentage of dark matter is %d.%d%%.\\n", int_part, dec_part1);

    return 0;
}
"""
    # The final prompt requires printing the calculated percentage and memory usage.
    # I am printing the C code for completeness as requested by the persona.
    print("--- Hypothetical Wuxing C Code ---")
    print(c_code)
    print("--- Analysis Results ---")
    print(f"p (Percentage of dark matter): {p}%")
    print(f"z (Estimated memory usage): {z} D")

    # The final required output format is p:z
    print("\n--- Final Answer for Submission ---")
    print(f"{p}:{z}")


generate_wuxing_c_code()