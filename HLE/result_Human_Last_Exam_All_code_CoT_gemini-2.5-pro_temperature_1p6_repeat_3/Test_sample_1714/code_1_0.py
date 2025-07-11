import math

def solve_pandora_problem():
    """
    This function performs the step-by-step calculation and prints the final result
    along with the hypothetical C code for the Wuxing architecture.
    """

    print("Step 1: Write the hypothetical C code for the Wuxing architecture.\n")
    
    # The hypothetical C code for the problem.
    c_code = """/*
 * =====================================================================================
 *
 *       Filename:  dark_matter.c
 *
 *    Description:  Calculates the dark matter percentage of Pandora galaxy
 *                  on the Wuxing decimal architecture.
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <frac.h> // Hypothetical header for the 'frac' type

int main() {
    // --- Variable Declaration ---
    // A 'frac' variable uses 3 chars (2D each) for n, d, and e.
    // Total memory for one frac = 3 * 2D = 6D.
    // We declare 7 frac variables, so total memory z = 7 * 6D = 42D.
    
    frac v, r, G, m_luminous, m_total, m_dark, percentage;

    // --- Initialization ---
    // v = 200 km/s -> 2 * 10^2
    v = 2e2;
    // r = 10 kpc -> 1 * 10^1
    r = 1e1;
    // G ≈ 4.3e-6 -> (43/10) * 10^-6, fits into {n=43, d=10, e=-6}
    G = 43/10e-6;
    
    // --- Calculation ---
    // 1. Luminous Mass = (M/L Ratio) * (Luminosity Multiplier)
    // M/L = 3, L = 2e9. Result is 6 * 10^9
    m_luminous = 3 * 2e9;

    // 2. Total Mass: M_total = (v^2 * r) / G
    m_total = v * v * r / G;

    // 3. Dark Matter Mass
    m_dark = m_total - m_luminous;

    // 4. Dark Matter Percentage
    percentage = m_dark * 100 / m_total;
    
    // --- Output ---
    // We assume the special %f specifier can print the result correctly.
    // The final calculation will be presented by the host script.
    printf("Dark Matter Percentage: %f\\n", percentage); 
    
    return 0;
}
"""
    print("--- Wuxing C Code ---")
    print(c_code)
    print("---------------------\n")


    print("Step 2: Perform the calculation based on the C code's logic.\n")
    # Constants and given values
    v = 200.0  # km/s
    r = 10.0   # kpc
    G = 4.3e-6 # (kpc/M_sun)*(km/s)^2, approximated as in the C code

    # Luminous Mass = (M/L ratio) * Luminosity = 3 * 2e9
    m_luminous = 6e9

    # Total Mass = v^2 * r / G
    m_total = (v**2 * r) / G

    # Dark Mass = Total Mass - Luminous Mass
    m_dark = m_total - m_luminous

    # Percentage = (Dark Mass / Total Mass) * 100
    dark_matter_percentage = (m_dark / m_total) * 100
    
    # Round percentage to 1 decimal place
    p_rounded = round(dark_matter_percentage, 1)

    # Memory usage calculation
    num_frac_variables = 7
    size_of_char_in_D = 2
    z_memory_usage = num_frac_variables * 3 * size_of_char_in_D

    print("--- Calculation Details ---")
    print(f"Total Mass Equation: M_total = (v^2 * r) / G")
    print(f"M_total = ({v:.0f}^2 * {r:.0f}) / {G:.1e} = {m_total:.4e} solar masses")
    print(f"\nLuminous Mass: M_lum = 3 * 2e9 = {m_luminous:.4e} solar masses")
    print(f"\nDark Mass: M_dark = M_total - M_lum = {m_dark:.4e} solar masses")
    print(f"\nPercentage Equation: p = (M_dark / M_total) * 100")
    print(f"p = ({m_dark:.4e} / {m_total:.4e}) * 100 = {dark_matter_percentage:.2f}%")
    print(f"\nRounded Percentage (p): {p_rounded}%")
    print(f"Memory Usage (z): {z_memory_usage}D")

    print("\nStep 3: Format the final answer as p:z.\n")
    final_answer = f"{p_rounded}:{z_memory_usage}"
    print(f"Final Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")


solve_pandora_problem()