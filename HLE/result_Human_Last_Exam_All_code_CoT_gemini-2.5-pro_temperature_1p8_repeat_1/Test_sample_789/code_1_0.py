import math

def main():
    """
    Performs the mass calculation according to Titan 5-bit architecture rules
    and computes the final absolute error.
    """
    
    # --- Step 1 & 2: Define initial fractions ---
    # Mass = (4/3) * pi * r^3 * density
    f_4_3 = (4, 3)
    f_pi = (22, 7)
    # r = 0.5 -> 1/2, so r^3 = 1/8
    f_r3 = (1, 8)
    # density = 0.9 -> 9/10
    f_rho = (9, 10)
    
    print("Titan Calculation Derivation for Mass of Rock")
    print("Formula: Mass = (4/3) * pi * r^3 * rho")
    print("Initial fraction substitution:")
    print(f"Mass = {f_4_3[0]}/{f_4_3[1]} * {f_pi[0]}/{f_pi[1]} * {f_r3[0]}/{f_r3[1]} * {f_rho[0]}/{f_rho[1]}\n")

    # --- Step 3: Re-group calculation to manage overflows ---
    print("Re-grouping to manage intermediate values:")
    print("Mass = ((4/3) * (9/10)) * ((22/7) * (1/8))\n")

    # --- Step 4: Calculate the first group ---
    # A = (4/3) * (9/10) -> (4*9)/(3*10) = 36/30 -> overflow
    # Simplify first: (4/10) -> 2/5 and (9/3) -> 3/1. (2/5)*(3/1) = 6/5
    group1_result = (6, 5)
    print("Calculating first group: (4/3) * (9/10)")
    print("Simplifying terms first gives (4*3)/10 = 12/10, which reduces to 6/5.")
    print(f"Result of first group: {group1_result[0]}/{group1_result[1]}\n")

    # --- Step 5: Calculate the second group and handle overflow ---
    # B = (22/7) * (1/8) = 22/56. Denominator > 31. OVERFLOW.
    group2_val = 22 / 56
    # We must approximate. 22/56 ~= 0.3928.
    # The fraction 7/18 (~0.3888) is chosen as it's a close approximation 
    # that allows the next calculation step to succeed.
    group2_approx = (7, 18)
    print("Calculating second group: (22/7) * (1/8)")
    print("This results in 22/56. The denominator 56 is > 31 (overflow).")
    print(f"We must replace 22/56 (value ~{group2_val:.4f}) with a valid 5-bit fraction.")
    print(f"The chosen approximation is {group2_approx[0]}/{group2_approx[1]} (value ~{group2_approx[0]/group2_approx[1]:.4f}).\n")

    # --- Step 6: Final calculation ---
    # Final = (6/5) * (7/18)
    # Simplify 6/18 -> 1/3.
    # (1/5) * (7/3) = 7/15
    final_result_fraction = (7, 15)
    print("Calculating final result: (Result of group 1) * (Approximation of group 2)")
    print(f"Mass = {group1_result[0]}/{group1_result[1]} * {group2_approx[0]}/{group2_approx[1]}")
    print("Simplifying before multiplying: 6/18 reduces to 1/3.")
    print(f"Mass = (1/5) * (7/3) = {final_result_fraction[0]}/{final_result_fraction[1]}\n")

    # --- Step 7: Error Calculation ---
    calculated_mass = final_result_fraction[0] / final_result_fraction[1]
    # True mass = (4/3) * pi * (0.5)^3 * 0.9 = 0.15 * pi
    true_mass = 0.15 * math.pi
    
    absolute_error = abs(true_mass - calculated_mass)
    # The question asks for the error rounded to 0.001
    rounded_error = round(absolute_error, 3)

    print("Error Calculation:")
    print(f"Final calculated mass: {calculated_mass:.6f}")
    print(f"True mass (0.15 * pi): {true_mass:.6f}")
    print(f"Absolute error: {absolute_error:.6f}")
    print(f"Smallest absolute error rounded to 0.001 (e): {rounded_error}")

if __name__ == "__main__":
    main()
