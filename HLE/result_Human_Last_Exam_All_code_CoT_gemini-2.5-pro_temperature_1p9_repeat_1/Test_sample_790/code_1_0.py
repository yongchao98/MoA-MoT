def check_overflow(n, term_name):
    """Checks if a number exceeds the 5-bit limit (0-31)."""
    if n > 31:
        # This operation would not be possible on Titan's hardware.
        # We stop this calculation path.
        raise ValueError(f"Constraint Violation: {term_name} {n} exceeds 31.")
    print(f"  - {term_name} {n} is valid.")

def print_frac(n, d, name):
    """Prints a fraction."""
    print(f"{name} = {n}/{d}")

def main():
    """
    Solves the physics problem using Titan's 5-bit fractional arithmetic.
    """
    print("--- Titan Superconducting Computer Simulation ---")
    print("Goal: Calculate the required force F.\n")
    print("Derived formula: F = 2 * g * m * sqrt(2)\n")

    # --- Step 1: Define Constants as 5-bit Fractions ---
    print("Step 1: Defining constants with valid fractional approximations.")
    # For pi, using 22/7 or other precise values leads to overflow during mass calculation.
    # We must use a cruder approximation.
    pi_n, pi_d = 3, 1
    print_frac(pi_n, pi_d, "Pi approximation (pi)")

    # For gravity g (9.8 m/s^2), using 29/3 (9.67) leads to overflow.
    # We must use g=10.
    g_n, g_d = 10, 1
    print_frac(g_n, g_d, "Gravity approximation (g)")

    # Physical parameters from the problem
    rho_n, rho_d = 9, 10  # 0.9 kg/cm^3
    print_frac(rho_n, rho_d, "Density (rho)")
    r_n, r_d = 1, 2      # 0.5 cm
    print_frac(r_n, r_d, "Radius (r)")
    const_4_3_n, const_4_3_d = 4, 3
    print("")

    # --- Step 2: Calculate Mass (m) = rho * (4/3) * pi * r^3 ---
    print("Step 2: Calculating mass (m) using fractional arithmetic.")
    # r^3 = (1/2)^3 = 1/8
    r3_n, r3_d = 1, 8

    # m_part1 = rho * (4/3) = (9/10) * (4/3)
    # Cross-cancellation: 9/3 = 3, 4/10 -> 2/5
    m1_n = (9 // 3) * (4 // 2)
    m1_d = (10 // 2) * (3 // 3)
    print("Intermediate m_part1 = (9/10)*(4/3) -> (3/1)*(2/5) = 6/5")

    # m_part2 = m_part1 * pi = (6/5) * (3/1)
    # Check intermediate product for overflow before simplifying
    check_overflow(6 * 3, "Numerator (6*3)")
    check_overflow(5 * 1, "Denominator (5*1)")
    m2_n, m2_d = 18, 5
    print("Intermediate m_part2 = (6/5)*(3/1) = 18/5")

    # m_final = m_part2 * r^3 = (18/5) * (1/8)
    # Check intermediate product for overflow
    check_overflow(18 * 1, "Numerator (18*1)")
    # 5*8=40 would overflow. We MUST use cross-cancellation.
    print("  - Denominator (5*8) = 40. This would overflow.")
    print("  - Applying algebraic simplification (cross-cancellation of 18 and 8 by 2):")
    m_final_n = (18 // 2) * 1
    m_final_d = 5 * (8 // 2)
    print(f"  - Simplified to ({m_final_n}/5) * (1/{m_final_d//5})")
    check_overflow(m_final_n, "Final mass numerator")
    check_overflow(m_final_d, "Final mass denominator")
    print_frac(m_final_n, m_final_d, "\nFinal mass (m)")
    print("")

    # --- Step 3: Calculate Force F = 2 * g * m * sqrt(2) ---
    print("Step 3: Calculating force (F).")
    # Part 1: (2/1) * g * m = (2/1) * (10/1) * (9/20)
    f_part1 = (2 * 10 * 9) // (1 * 1 * 20) # Simplified value
    print("Intermediate F_part = (2/1)*(10/1)*(9/20) simplifies to 9/1")

    # Part 2: F = F_part * sqrt(2) = 9 * sqrt(2)
    print("\nCalculation is now F = 9 * sqrt(2)")
    # Using a precise sqrt(2) like 24/17 results in 9*24 = 216 -> overflow.
    # We must use expansion and approximation. Let sqrt(2) approx 7/5 = 1 + 2/5
    sqrt2_approx_n, sqrt2_approx_d = 7, 5
    print_frac(sqrt2_approx_n, sqrt2_approx_d, "Approximation for sqrt(2)")
    
    # F = 9 * (7/5) = 9 * (1 + 2/5) = 9 + 9 * (2/5) = 9 + 18/5
    term2_n, term2_d = 18, 5
    print(f"Expansion leads to F = 9/1 + {term2_n}/{term2_d}")

    # Check if adding them directly is possible: (9*5 + 18*1)/(1*5) = (45+18)/5 = 63/5.
    print("Directly adding 9/1 and 18/5 causes overflow (45+18=63).")
    
    # Apply Rule 5: Simplify by using a less precise approximation for the smaller term.
    # 18/5 = 3.6. We approximate this with the closest integer, 4.
    term2_approx_n, term2_approx_d = 4, 1
    print(f"Approximating the smaller term {term2_n}/{term2_d} (3.6) with {term2_approx_n}/{term2_approx_d} (4.0).")
    
    # Final addition: F = 9/1 + 4/1
    F_final_n = 9 + 4
    F_final_d = 1
    print("\n--- Final Calculation ---")
    print(f"F = 9 / 1 + {term2_approx_n} / {term2_approx_d}")
    print("The final equation is:")
    print(f"Force = {F_final_n} / {F_final_d}")

if __name__ == '__main__':
    main()
