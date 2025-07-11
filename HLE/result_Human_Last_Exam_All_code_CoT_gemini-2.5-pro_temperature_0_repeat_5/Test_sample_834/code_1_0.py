import math

def solve_beam_problem():
    """
    This script solves the beam deflection problem by following these steps:
    1. Derives the formula for the force F required to make the tip deflection zero using the principle of superposition.
    2. Calculates the moments of inertia I_zz and I_ss for the given cross-section.
    3. Uses the provided formulas and calculated moments of inertia to find the numerical values for L and q0.
    4. Substitutes these values into the formula for F to get the final result.
    """

    # Step 1: Define geometric and material properties
    # a = 12^(1/4)
    a_pow_4 = 12

    # Step 2: Calculate moments of inertia I_zz and I_ss
    # For a large square (3a x 3a) with two smaller squares (a x a) removed.
    # I_zz = I_zz_main - 2 * (I_zz_cutout_centroid + A_cutout * d_s^2)
    # I_zz_main = (3a)*(3a)^3/12 = 27*a^4/4
    # I_zz_cutout = a^4/12 + a^2*(a/2)^2 = a^4/3
    # I_zz = 27*a^4/4 - 2*(a^4/3) = (81-8)*a^4/12 = 73*a^4/12
    I_zz = 73 * a_pow_4 / 12

    # I_ss = I_ss_main - 2 * (I_ss_cutout_centroid + A_cutout * d_z^2)
    # I_ss_main = (3a)*(3a)^3/12 = 27*a^4/4
    # I_ss_cutout = a^4/12 + a^2*(a)^2 = 13*a^4/12
    # I_ss = 27*a^4/4 - 2*(13*a^4/12) = (81-26)*a^4/12 = 55*a^4/12
    I_ss = 55 * a_pow_4 / 12

    # Step 3: Calculate L and q0 from the given data
    L = (30 * I_zz) / 73
    q0 = (9 * I_ss) / 55

    # Step 4: Calculate the required force F
    # The derived formula for F is F = (13 * q0 * L) / 270
    F = (13 * q0 * L) / 270

    # Step 5: Print the results
    print("--- Intermediate Calculations ---")
    print(f"Value of a^4: {a_pow_4}")
    print(f"Moment of inertia I_zz: {I_zz}")
    print(f"Moment of inertia I_ss: {I_ss}")
    print(f"Length L: {L}")
    print(f"Maximum distributed load q0: {q0}")
    print("\n--- Final Calculation ---")
    print("The force F is calculated using the formula: F = (13 * q0 * L) / 270")
    # The problem asks to output each number in the final equation
    print(f"Substituting the values, the final equation is: F = (13 * {q0} * {L}) / 270")
    print(f"The required force F is: {F}")
    
    # Final answer in the required format
    print(f"\n<<<{int(round(F))}>>>")

solve_beam_problem()