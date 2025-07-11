import math

def solve_star_angle_problem():
    """
    Solves the relativistic star angle problem by calculating the ratio of Doppler factors.
    """

    # --- Step 1: Set up the equations using data from Frame 2 ---
    print("Step 1: Set up equations for Doppler factors based on Frame 2 data.\n")
    # In Frame 1, for any two stars i,j, cos(theta_ij) = -1/3.
    # The Lorentz-invariant relation gives: 1 - cos(theta'_ij) = (4/3) / (D_i * D_j)
    # We are given the following angles in Frame 2:
    # Angle between S1 and S2: pi/2 -> cos(theta'_12) = 0
    # Angle between S1 and S3: 3*pi/4 -> cos(theta'_13) = -1/sqrt(2)
    # Angle between S2 and S3: 3*pi/4 -> cos(theta'_23) = -1/sqrt(2)
    
    cos_theta_12_prime = 0
    cos_theta_13_prime = -1 / math.sqrt(2)

    # Equation for D1 * D2
    # 1 - cos(theta'_12) = (4/3) / (D1 * D2) => D1 * D2 = (4/3) / (1 - 0)
    val_D1_D2 = 4/3
    print(f"From angle theta'_12, we get the equation: D_1 * D_2 = 4/3 = {val_D1_D2:.6f}")

    # Equation for D1 * D3
    # 1 - cos(theta'_13) = (4/3) / (D1 * D3) => D1 * D3 = (4/3) / (1 - (-1/sqrt(2)))
    val_1_minus_cos_13 = 1 - cos_theta_13_prime
    val_D1_D3 = (4/3) / val_1_minus_cos_13
    print(f"From angle theta'_13, we get the equation: D_1 * D_3 = (4/3) / (1 - cos(3pi/4)) = {val_D1_D3:.6f}")

    # Equation for D2 * D3
    val_D2_D3 = val_D1_D3  # Since theta'_23 is the same as theta'_13
    print(f"From angle theta'_23, we get the equation: D_2 * D_3 = (4/3) / (1 - cos(3pi/4)) = {val_D2_D3:.6f}\n")

    # --- Step 2: Solve the system of equations for the ratio D3/D1 ---
    print("Step 2: Solve the system for the ratio D_3 / D_1.\n")
    # From D_1*D_3 = D_2*D_3, we can conclude that D_1 = D_2.
    print("Since D_1*D_3 = D_2*D_3, we know that D_1 = D_2.")
    
    # Substitute D_1 = D_2 into the first equation: D_1^2 = 4/3
    # D1 = sqrt(4/3) = 2/sqrt(3)
    D1 = math.sqrt(val_D1_D2)
    print(f"Substituting into the first equation: D_1^2 = 4/3, so D_1 = sqrt(4/3) = {D1:.6f}")
    
    # Now we find D_3 using D_1*D_3
    # D3 = (val_D1_D3) / D1
    D3 = val_D1_D3 / D1
    print(f"Using the second equation: D_3 = (D_1*D_3) / D_1 = {val_D1_D3:.6f} / {D1:.6f} = {D3:.6f}\n")
    
    # --- Step 3: Calculate the final ratio ---
    print("Step 3: Calculate the final requested value.\n")
    # The value to find is R = (1 - cos(theta'_14)) / (1 - cos(theta'_34)), which simplifies to D_3 / D_1.
    R = D3 / D1
    
    # Print the numbers in the final calculation
    print("The final equation is Ratio = D_3 / D_1.")
    print(f"The number for the numerator is D_3 = {D3:.6f}")
    print(f"The number for the denominator is D_1 = {D1:.6f}")
    
    # --- Final Answer ---
    final_numeric_value = 2 - math.sqrt(2)
    print("\nResult:")
    print(f"The ratio (1 - cos(theta_14)) / (1 - cos(theta_34)) is equal to D_3 / D_1 = {R:.6f}")
    print(f"The exact symbolic value of the ratio is 2 - sqrt(2).")
    print(f"The final numeric value is: {final_numeric_value}")


solve_star_angle_problem()