import math

def solve_titan_mass_problem():
    """
    This function demonstrates the derivation of the rock's mass calculation
    under the constraints of the Titan 5-bit architecture.
    """

    # Helper function to check if a number fits in 5 bits (0-31)
    def is_5_bit(n):
        return 0 <= n <= 31

    # Helper function to simplify a fraction n/d
    def simplify(n, d):
        if d == 0:
            return float('inf'), float('inf')
        common = math.gcd(n, d)
        return n // common, d // common

    print("Titan-V Fractional Calculation for Mass of a Rock")
    print("="*50)

    # Step 1: Define initial values as 5-bit fractions.
    # A 5-bit integer can hold values from 0 to 2^5 - 1 = 31.
    print("1. Representing all values as 5-bit fractions:")
    r_num, r_den = 1, 2       # r = 0.5 cm
    rho_num, rho_den = 9, 10    # rho = 0.9 kg/cm^3
    pi_num, pi_den = 22, 7      # pi is approximated as 22/7
    c_num, c_den = 4, 3       # constant 4/3

    print(f"   Radius (r) = 0.5 cm              -> {r_num}/{r_den}")
    print(f"   Density (rho) = 0.9 kg/cm^3      -> {rho_num}/{rho_den}")
    print(f"   Pi (pi)                          -> {pi_num}/{pi_den} (approximation)")
    print(f"   Constant                         -> {c_num}/{c_den}")
    
    # Calculate r^3
    r3_num, r3_den = r_num**3, r_den**3
    
    print("\n2. Assembling the equation for mass:")
    print(f"   Mass = ({c_num}/{c_den}) * ({pi_num}/{pi_den}) * ({r_num}/{r_den})^3 * ({rho_num}/{rho_den})")
    print(f"   Mass = ({c_num}/{c_den}) * ({pi_num}/{pi_den}) * ({r3_num}/{r3_den}) * ({rho_num}/{rho_den})")

    print("\n3. Performing calculation, simplifying at each step to avoid overflow.")
    # We strategically reorder to (4/3)*(9/10) * (1/8) * (22/7)
    
    # Term 1: (4/3) * (9/10)
    print("\n   Step A: (4/3) * (9/10)")
    n1_unsimplified, d1_unsimplified = c_num * rho_num, c_den * rho_den
    print(f"      Directly multiplying gives {n1_unsimplified}/{d1_unsimplified} = 36/30. The numerator 36 overflows the 5-bit limit.")
    n1, d1 = simplify(n1_unsimplified, d1_unsimplified)
    print(f"      Simplifying gives a valid fraction: {n1}/{d1}")

    # Term 2: (6/5) * (1/8)
    print(f"\n   Step B: ({n1}/{d1}) * (1/8)")
    n2_unsimplified, d2_unsimplified = n1 * r3_num, d1 * r3_den
    print(f"      Directly multiplying gives {n2_unsimplified}/{d2_unsimplified} = 6/40. The denominator 40 overflows the 5-bit limit.")
    n2, d2 = simplify(n2_unsimplified, d2_unsimplified)
    print(f"      Simplifying gives a valid fraction: {n2}/{d2}")

    # Term 3: (3/20) * (22/7)
    print(f"\n   Step C: ({n2}/{d2}) * (22/7)")
    n3_unsimplified, d3_unsimplified = n2 * pi_num, d2 * pi_den
    print(f"      Multiplying gives {n3_unsimplified}/{d3_unsimplified} = 66/140. Both numbers overflow the 5-bit limit.")
    final_num, final_den = simplify(n3_unsimplified, d3_unsimplified)
    print(f"      Simplifying gives {final_num}/{final_den}. This is the most accurate result, but still overflows (33 > 31).")

    print("\n4. Finding the best 5-bit approximation for the invalid result.")
    target_val = final_num / final_den
    print(f"   The target value is {final_num}/{final_den} ≈ {target_val:.5f}.")

    best_n, best_d = -1, -1
    min_error_approx = float('inf')

    for b in range(1, 32):
        a = round(target_val * b)
        if is_5_bit(a) and a > 0:
            current_error = abs(target_val - a / b)
            if current_error < min_error_approx:
                min_error_approx = current_error
                best_n, best_d = a, b
    
    best_n, best_d = simplify(best_n, best_d)

    print(f"   Searching all fractions a/b (a,b <= 31), the closest one is {best_n}/{best_d}.")
    
    print("\n5. Final Calculation Summary:")
    print("   The full equation is: (4/3) * (22/7) * (1/8) * (9/10)")
    print(f"   This calculates to 33/70, which is not representable on Titan.")
    print(f"   It must be approximated by the closest valid fraction, which is {best_n}/{best_d}.")
    print(f"   The final calculated mass is therefore {best_n}/{best_d}.")

    # Final step is to calculate the smallest possible absolute error.
    true_mass = (4/3) * math.pi * (0.5**3) * 0.9
    our_mass = best_n / best_d
    absolute_error = abs(true_mass - our_mass)
    
    print("\n6. Calculating Smallest Absolute Error:")
    print(f"   True Mass (using high-precision pi) ≈ {true_mass:.6f}")
    print(f"   Our Calculated Mass ({best_n}/{best_d})         ≈ {our_mass:.6f}")
    print(f"   Absolute Error                        = |{true_mass:.6f} - {our_mass:.6f}| ≈ {absolute_error:.6f}")

    # The final answer is the error rounded to 0.001 (3 decimal places)
    final_answer = round(absolute_error, 3)
    print(f"\n   The error, rounded to the nearest 0.001, is: {final_answer}")


if __name__ == '__main__':
    solve_titan_mass_problem()
