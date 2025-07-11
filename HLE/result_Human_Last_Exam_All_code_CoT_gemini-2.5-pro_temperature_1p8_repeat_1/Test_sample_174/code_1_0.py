import math
from fractions import Fraction

def calculate_orbifold_euler_characteristic():
    """
    Calculates the orbifold Euler characteristic of the quotient stack of smooth plane quartics.
    """
    # Genus of a smooth plane quartic curve
    g = 3

    # --- Step 1: Calculate chi_orb(M_3) ---
    # The formula is chi_orb(M_g) = B_{2g} / (2g * (2g - 2)).
    # We need B_6 = 1/42.
    B_6 = Fraction(1, 42)
    chi_M_3_denominator = 2 * g * (2 * g - 2)
    chi_M_3 = B_6 / chi_M_3_denominator

    # --- Step 2: Calculate chi_orb(H_3) ---
    # First, calculate chi_orb(M_{0, 2g+2}) for g=3, which is chi_orb(M_{0, 8}).
    n_points = 2 * g + 2  # This is 8
    # Formula: chi_orb(M_{0,n}) = (-1)^(n-3) * (n-3)!
    n_minus_3 = n_points - 3
    chi_M_0_8 = (-1)**n_minus_3 * math.factorial(n_minus_3)
    
    # Then, calculate chi_orb(H_3)
    # Formula: chi_orb(H_g) = (1/2) * chi_orb(M_{0, 2g+2})
    chi_H_3 = Fraction(1, 2) * chi_M_0_8

    # --- Step 3: Final Calculation ---
    # chi_orb([U/G]) = chi_orb(M_3) - chi_orb(H_3)
    result = chi_M_3 - chi_H_3
    
    print("The calculation proceeds as follows:")
    print("chi_orb([U/G]) = chi_orb(M_3^nh) = chi_orb(M_3) - chi_orb(H_3)")
    print("\nCalculation of chi_orb(M_3):")
    print(f"chi_orb(M_3) = B_6 / (2*3 * (2*3 - 2)) = (1/42) / {chi_M_3_denominator} = {chi_M_3}")
    
    print("\nCalculation of chi_orb(H_3):")
    print(f"chi_orb(M_{0, 8}) = (-1)^({n_points}-3) * ({n_points}-3)! = (-1)^{n_minus_3} * {n_minus_3}! = {chi_M_0_8}")
    print(f"chi_orb(H_3) = (1/2) * chi_orb(M_{0, 8}) = (1/2) * {chi_M_0_8} = {chi_H_3}")

    print("\nFinal Result:")
    # To satisfy the output format requirement, print each number in the final equation.
    print(f"chi_orb([U/G]) = {chi_M_3} - ({chi_H_3}) = {result}")

    # Output the final fraction clearly
    final_numerator = result.numerator
    final_denominator = result.denominator
    print(f"\nThe orbifold Euler characteristic of [U/G] is {final_numerator}/{final_denominator}.")


calculate_orbifold_euler_characteristic()
<<<60481/1008>>>