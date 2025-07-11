import math

def solve_star_angle_problem():
    """
    Solves the relativistic star angle problem by applying Lorentz invariants.

    The solution proceeds in the following steps:
    1.  Determine the angular separation of the stars in the first reference frame.
    2.  Use the Lorentz invariant quantity related to the 4-momenta of photons from the stars.
    3.  Set up a system of equations based on the observer data in the second frame.
    4.  Solve for the ratio of photon energies in the second frame.
    5.  Relate this energy ratio to the desired ratio of angle cosines.
    6.  Calculate and print the final result.
    """
    print("Step 1: Analyzing the first reference frame")
    # In the first frame, the four stars S1, S2, S3, S4 form a regular tetrahedron
    # on the celestial sphere. Let k_i be the unit vector pointing to star S_i.
    # Due to the tetrahedral symmetry, k_1 + k_2 + k_3 + k_4 = 0.
    # The angle 'alpha' between any two stars is the same. cos(alpha) = k_i . k_j for i != j.
    # Dotting the vector sum with k_1 gives: k_1.(k_1 + k_2 + k_3 + k_4) = 0
    # k_1.k_1 + k_1.k_2 + k_1.k_3 + k_1.k_4 = 0  => 1 + 3*cos(alpha) = 0
    cos_alpha = -1/3
    print(f"In the first frame, the angle 'alpha' between any pair of stars is constant.")
    print(f"From the tetrahedral geometry, we find cos(alpha) = -1/3 = {cos_alpha:.4f}.")
    one_minus_cos_alpha = 1 - cos_alpha
    print(f"Thus, the term (1 - cos(theta_ij)) is 1 - (-1/3) = {one_minus_cos_alpha:.4f} for all pairs.\n")

    print("Step 2: Using the Lorentz Invariant")
    # The scalar product of the 4-momenta of two photons (p_i, p_j) is a Lorentz invariant.
    # p_i . p_j = (E_i*E_j / c^2) * (1 - k_i . k_j)
    # This means (E_i * E_j * (1 - cos(theta_ij))) is the same in all inertial frames.
    # Since E_i and (1 - cos(theta_ij)) are the same for all pairs in Frame 1, the product
    # E'_i * E'_j * (1 - cos(theta'_ij)) must be a constant in Frame 2.
    print("The quantity E'_i * E'_j * (1 - cos(theta'_ij)) is a constant for any pair of stars (i,j) in the second frame.\n")

    print("Step 3: Setting up equations for the second frame")
    # Let E'_1, E'_2, E'_3, E'_4 be the observed energies in the second frame.
    # Given data for Frame 2:
    # theta'_12 = pi/2         => cos(theta'_12) = 0
    # theta'_13 = 3*pi/4       => cos(theta'_13) = -1/sqrt(2)
    # theta'_23 = 3*pi/4       => cos(theta'_23) = -1/sqrt(2)
    one_minus_cos_12 = 1.0
    one_minus_cos_13 = 1.0 - (-1.0 / math.sqrt(2))
    
    print("Based on the given angles, we have:")
    print(f"1 - cos(theta'_12) = {one_minus_cos_12:.4f}")
    print(f"1 - cos(theta'_13) = 1 - cos(theta'_23) = {one_minus_cos_13:.4f}")
    print("\nThis gives the following system of equations (where C is a constant):")
    print(f"(1) E'_1 * E'_2 * ({one_minus_cos_12:.4f}) = C")
    print(f"(2) E'_1 * E'_3 * ({one_minus_cos_13:.4f}) = C")
    print(f"(3) E'_2 * E'_3 * ({one_minus_cos_13:.4f}) = C\n")
    
    print("Step 4: Solving for the ratio of energies")
    print("From equations (2) and (3), by dividing them, we deduce that E'_1 = E'_2.")
    print("Substituting E'_1 = E'_2 into (1), we get: (E'_1)^2 = C.")
    # Substituting C = (E'_1)^2 into (2): E'_1 * E'_3 * (1 + 1/sqrt(2)) = (E'_1)^2
    # This gives the ratio E'_3 / E'_1
    ratio_E3_E1 = 1.0 / one_minus_cos_13
    print(f"Substituting C into (2) and solving for E'_3 / E'_1 gives: 1 / (1 + 1/sqrt(2)) = {ratio_E3_E1:.4f}\n")

    print("Step 5: Finding the final required ratio")
    # We need to find the value of (1 - cos(theta'_14)) / (1 - cos(theta'_34))
    # Using the invariance principle:
    # 1 - cos(theta'_14) = C / (E'_1 * E'_4)
    # 1 - cos(theta'_34) = C / (E'_3 * E'_4)
    # The ratio is: (C / (E'_1 * E'_4)) / (C / (E'_3 * E'_4)) = E'_3 / E'_1
    print("The required ratio (1 - cos(theta'_14)) / (1 - cos(theta'_34)) simplifies exactly to the energy ratio E'_3 / E'_1.\n")
    
    print("Step 6: Final Calculation")
    # Let's simplify the expression 1 / (1 + 1/sqrt(2))
    # 1 / (1 + 1/sqrt(2)) = 1 / ((sqrt(2)+1)/sqrt(2)) = sqrt(2)/(sqrt(2)+1)
    # Rationalize the denominator by multiplying by (sqrt(2)-1)/(sqrt(2)-1):
    # (sqrt(2) * (sqrt(2)-1)) / ((sqrt(2)+1)*(sqrt(2)-1)) = (2 - sqrt(2)) / (2 - 1) = 2 - sqrt(2)
    final_value_symbolic = "2 - sqrt(2)"
    final_value_numeric = 2.0 - math.sqrt(2)
    print("The final result is the value of the ratio we found.")
    print(f"The final equation is: (1 - cos(theta_14)) / (1 - cos(theta_34)) = {final_value_symbolic}")
    print(f"The numbers in this final equation are {2} and sqrt({2}).")
    print(f"The numerical value of the expression is: {final_value_numeric}")

solve_star_angle_problem()