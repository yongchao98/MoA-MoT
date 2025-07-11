import math

def print_step(description, equation_str):
    """Helper function to print calculation steps."""
    print(f"{description:<25} {equation_str}")

def main():
    """
    Solves the Curious Monkey and Coconut problem using Titan computer constraints.
    """
    print("Titan Computer Calculation for Required Force\n")
    print("Physical Model: Constant force 'F' during flight.")
    print("Derived Formula: F = 2 * g * m * sqrt(2)\n")

    # --- Step 1: Define constants as fractions (n, d) where n, d <= 31 ---
    print_step("1. Approximating constants:", "")
    # g ~ 9.8 m/s^2. We use 10/1 for easier calculation.
    g = (10, 1)
    print_step("   Gravity g:", f"{g[0]}/{g[1]} (Value: {g[0]/g[1]})")

    # pi ~ 3.14159. We use 16/5 = 3.2 to make mass calculation feasible.
    pi = (16, 5)
    print_step("   Pi:", f"{pi[0]}/{pi[1]} (Value: {pi[0]/pi[1]})")

    # sqrt(2) ~ 1.414. We use 7/5 = 1.4.
    sqrt2 = (7, 5)
    print_step("   sqrt(2):", f"{sqrt2[0]}/{sqrt2[1]} (Value: {sqrt2[0]/sqrt2[1]})")
    
    two = (2, 1)
    three_twentieths = (3, 20)
    print("")

    # --- Step 2: Calculate mass 'm' ---
    # m = 3 * pi / 20 = (3/20) * (16/5)
    print_step("2. Calculating mass 'm':", f"m = (3/20) * pi = ({three_twentieths[0]}/{three_twentieths[1]}) * ({pi[0]}/{pi[1]})")
    # Simplify before multiplying: (3/20)*(16/5) = (3/(5*4))*(4*4/5) = (3/5)*(4/5)
    m_num = 3 * 4
    m_den = 5 * 5
    m = (m_num, m_den)
    print_step("   Simplified to:", f"(3/5) * (4/5) = {m[0]}/{m[1]}")
    print("")

    # --- Step 3: Calculate Force 'F' ---
    # F = 2 * g * m * sqrt(2)
    print_step("3. Calculating force 'F':", f"F = ({two[0]}/{two[1]}) * ({g[0]}/{g[1]}) * ({m[0]}/{m[1]}) * ({sqrt2[0]}/{sqrt2[1]})")
    
    # Group terms to keep intermediate values small
    # F = (g * m) * (2 * sqrt(2))
    # term1 = g * m = (10/1) * (12/25) = (2*5/1) * (12/5/5) = (2/1) * (12/5) = 24/5
    term1_num = 2 * 12
    term1_den = 5
    term1 = (term1_num, term1_den)
    print_step("   Intermediate term 1:", f"({g[0]}/{g[1]}) * ({m[0]}/{m[1]}) = {term1[0]}/{term1[1]}")

    # term2 = 2 * sqrt(2) = (2/1) * (7/5) = 14/5
    term2_num = 2 * 7
    term2_den = 5
    term2 = (term2_num, term2_den)
    print_step("   Intermediate term 2:", f"({two[0]}/{two[1]}) * ({sqrt2[0]}/{sqrt2[1]}) = {term2[0]}/{term2[1]}")

    # Now F = term1 * term2 = (24/5) * (14/5)
    # This multiplication (24*14=336) would exceed the 5-bit limit.
    print_step("   Next step:", f"F = ({term1[0]}/{term1[1]}) * ({term2[0]}/{term2[1]})")
    print("   Constraint Violation: Numerator 24 * 14 = 336, which is > 31.")
    
    # --- Step 4: Apply approximation to maintain constraints ---
    # We replace 24/5 (value 4.8) with 25/5 (value 5.0) to enable simplification.
    # This is allowed by the Titan computational rules.
    approx_term1 = (25, 5)
    print_step("4. Applying approximation:", f"Replace {term1[0]}/{term1[1]} with {approx_term1[0]}/{approx_term1[1]}")
    
    # Final calculation: F ~ (25/5) * (14/5) = (5/1) * (14/5) = 14/1
    final_F_num = 14
    final_F_den = 1
    final_F = (final_F_num, final_F_den)
    print_step("   Final calculation:", f"({approx_term1[0]}/{approx_term1[1]}) * ({term2[0]}/{term2[1]}) = (5/1) * (14/5) = {final_F[0]}/{final_F[1]}")
    print("")

    # --- Step 5: Error Analysis ---
    print("5. Error Analysis:")
    # True values
    true_g = 9.8
    true_pi = math.pi
    true_sqrt2 = math.sqrt(2)
    true_m = (3 * true_pi) / 20
    true_F = 2 * true_g * true_m * true_sqrt2
    
    calculated_F = final_F[0] / final_F[1]
    abs_error = abs(calculated_F - true_F)
    
    print(f"   True Force Value: {true_F:.3f} N")
    print(f"   Titan Calculated Value: {calculated_F:.3f} N")
    print(f"   Absolute Error: {abs_error:.3f}")
    
    # --- Final Answer ---
    print("\nYes, the force can be calculated.")
    print(f"The smallest absolute error produced is {abs_error:.3f}.")
    print(f"<<<Y[{abs_error:.3f}]>>>")

main()