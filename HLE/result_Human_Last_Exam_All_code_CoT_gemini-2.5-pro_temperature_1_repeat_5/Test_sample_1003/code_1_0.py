import math

def solve_star_angle_problem():
    """
    Solves the relativistic star angle problem by calculating the ratio
    (1 - cos(theta'_14)) / (1 - cos(theta'_34)).
    """

    # --- Step 1: Theoretical foundation ---
    # The key physical principle is that the quantity (omega_i * omega_j * (1 - cos(theta_ij)))
    # is a Lorentz invariant. Using the Doppler shift formula omega'_i = D_i * omega_i,
    # we derive the angle transformation law:
    # 1 - cos(theta'_ij) = (1 / (Di * Dj)) * (1 - cos(theta_ij))
    print("The problem is solved using the relativistic aberration of light.")
    print("The key relation is: 1 - cos(theta'_ij) = (1 / (Di * Dj)) * (1 - cos(theta_ij))\n")

    # --- Step 2: Information from Frame 1 ---
    # In Frame 1, the stars form a regular tetrahedron. The angle between any two
    # stars as seen from the center has cos(theta) = -1/3.
    cos_theta_frame1 = -1/3
    one_minus_cos_theta_frame1 = 1 - cos_theta_frame1
    print(f"In Frame 1, for any pair of stars, 1 - cos(theta_ij) = 1 - (-1/3) = {one_minus_cos_theta_frame1:.4f}")

    # --- Step 3: Express the target ratio in terms of Doppler factors Di ---
    # Target = (1 - cos(theta'_14)) / (1 - cos(theta'_34))
    #        = [ (1/(D1*D4)) * (4/3) ] / [ (1/(D3*D4)) * (4/3) ]
    #        = D3 / D1
    print("The expression to find, (1 - cos(theta'_14)) / (1 - cos(theta'_34)), simplifies to the ratio D3 / D1.\n")

    # --- Step 4: Use information from Frame 2 ---
    print("Using the angles given for Frame 2 to set up equations for the Doppler factors:")
    # From theta'_12 = pi/2, we have 1 - cos(theta'_12) = 1.
    # 1 = (1 / (D1*D2)) * (4/3)  => D1 * D2 = 4/3
    d1_d2 = 4/3
    print(f"From theta'_12 = 90 deg, we find D1*D2 = 4/3 = {d1_d2:.4f}")

    # From theta'_13 = 3pi/4, we have 1 - cos(theta'_13) = 1 + 1/sqrt(2).
    one_minus_cos_theta_13 = 1 + 1/math.sqrt(2)
    # (1 + 1/sqrt(2)) = (1 / (D1*D3)) * (4/3) => D1*D3 = (4/3)/(1 + 1/sqrt(2))
    d1_d3 = d1_d2 / one_minus_cos_theta_13
    print(f"From theta'_13 = 135 deg, we find D1*D3 = (4/3) / (1 + 1/sqrt(2)) = {d1_d3:.4f}")

    # From theta'_23 = 3pi/4, we get an identical equation for D2*D3, which implies D1=D2.
    print("From theta'_23 = 135 deg, we find D2*D3 is the same, which implies D1 = D2.\n")

    # --- Step 5: Calculate the final ratio D3/D1 ---
    print("To find the final ratio D3/D1, we divide the equation for D1*D3 by the equation for D1*D1 (which is D1*D2):")
    # Ratio = D3/D1 = (D1*D3) / (D1*D1) = d1_d3 / d1_d2
    final_ratio = d1_d3 / d1_d2
    print(f"Ratio = [ (4/3)/(1 + 1/sqrt(2)) ] / [ 4/3 ] = 1 / (1 + 1/sqrt(2))")
    
    # --- Step 6: Final Result ---
    # The exact value is 2 - sqrt(2)
    exact_value_str = "2 - sqrt(2)"
    numerical_value = 2 - math.sqrt(2)
    
    print(f"\nSimplifying this expression gives the exact answer: {exact_value_str}")
    print(f"The numerical value is approximately {numerical_value:.8f}\n")
    
    # Per instructions, outputting each number in the final equation "Value = 2 - sqrt(2)"
    print("The final equation is of the form: A - sqrt(B)")
    print(f"The number A is: 2")
    print(f"The number B under the square root is: 2")

solve_star_angle_problem()