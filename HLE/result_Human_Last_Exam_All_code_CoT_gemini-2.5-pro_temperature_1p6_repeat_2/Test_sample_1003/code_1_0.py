import math

def solve_relativistic_aberration():
    """
    Solves the problem by applying the principles of special relativity to the
    aberration of light from four distant stars.
    """
    print("### Step-by-Step Solution ###")

    # Step 1: Analyze the first reference frame (S)
    print("\nStep 1: Analyze the geometry in the first reference frame.")
    print("In frame S, the apparent angle between any pair of the four stars is equal.")
    print("This means the direction vectors to the stars form a regular tetrahedron.")
    print("For this configuration, the cosine of the angle 'theta' between any two direction vectors is -1/3.")
    cos_theta_S = -1/3
    one_minus_cos_theta_S = 1 - cos_theta_S
    print(f"Therefore, for any pair of stars i and j, 1 - cos(theta_ij) = 1 - ({cos_theta_S}) = {one_minus_cos_theta_S:.3f}")

    # Step 2: Formulate the relationship between the two frames
    print("\nStep 2: Relate the angles in the two frames.")
    print("From the Lorentz invariance of the 4-momentum scalar product, we have the relation:")
    print("1 - cos(theta_ij) = A_i * A_j * (1 - cos(theta'_ij))")
    print("where theta'_ij are the angles in the second frame, and A_k are factors related to the relative motion.")

    # Step 3: Simplify the desired ratio
    print("\nStep 3: Simplify the expression we need to calculate.")
    print("We need to find the value of R = (1 - cos(theta'_14)) / (1 - cos(theta'_34)).")
    print("Using the relation from Step 2:")
    print("R = [ (1 - cos(theta_14))/(A1*A4) ] / [ (1 - cos(theta_34))/(A3*A4) ]")
    print(f"Since (1 - cos(theta_14)) and (1 - cos(theta_34)) both equal {one_minus_cos_theta_S:.3f}, the expression simplifies to:")
    print("R = A3 / A1")

    # Step 4: Use the data from the second frame to find the ratio R
    print("\nStep 4: Use the given angles in the second frame to find R.")
    cos_theta_prime_12 = 0
    cos_theta_prime_13 = -1 / math.sqrt(2)
    print(f"Given data for the second frame (S'):")
    print(f"  - Angle between S1 and S2 is 90 degrees, so cos(theta'_12) = {cos_theta_prime_12}")
    print(f"  - Angle between S1 and S3 is 135 degrees (3pi/4), so cos(theta'_13) = {cos_theta_prime_13:.3f}")

    print("\nLet's write down the equations for pairs (1,2) and (1,3):")
    print(f"(a) For (1,2): {one_minus_cos_theta_S:.3f} = A1 * A2 * (1 - {cos_theta_prime_12})  => A1 * A2 = 4/3")
    print(f"(b) For (1,3): {one_minus_cos_theta_S:.3f} = A1 * A3 * (1 - ({cos_theta_prime_13:.3f}))")

    print("\nFrom the problem's symmetry, we can deduce A1=A2, which means from (a), A1^2 = 4/3.")
    print("Dividing equation (a) by (b) is complex. A simpler way is to equate the left-hand sides:")
    print("A1 * A2 = A1 * A3 * (1 - cos(theta'_13))")
    print("Substituting A2 = A1:")
    print("A1^2 = A1 * A3 * (1 - cos(theta'_13))")
    print("Dividing by A1 (which is non-zero) gives: A1 = A3 * (1 - cos(theta'_13))")
    print("So, the ratio R = A3 / A1 is: R = 1 / (1 - cos(theta'_13))")

    # Step 5: Calculate the final value
    print("\nStep 5: Calculate the final value.")
    num_1 = 1
    den_1 = 1
    den_2 = -1
    sqrt_val = 2
    cos_theta_p_13_val = den_2 / math.sqrt(sqrt_val)
    print(f"Substituting the value of cos(theta'_13) = -1/sqrt(2):")
    print(f"R = {num_1} / ({den_1} - ({cos_theta_p_13_val:.3f}))")
    print(f"Final Equation: R = {num_1} / ({den_1} + {abs(den_2)}/sqrt({sqrt_val}))")

    print("\nSimplifying the expression:")
    print("R = 1 / (1 + 1/sqrt(2)) = sqrt(2) / (sqrt(2) + 1)")
    print("Multiplying the numerator and denominator by (sqrt(2) - 1):")
    print("R = (sqrt(2) * (sqrt(2) - 1)) / ((sqrt(2) + 1) * (sqrt(2) - 1))")
    print("R = (2 - sqrt(2)) / (2 - 1)")
    
    term_1 = 2
    sqrt_val_2 = 2
    print(f"Simplified Final Equation: R = {term_1} - sqrt({sqrt_val_2})")
    
    final_value = 2 - math.sqrt(2)
    print(f"\nThe numerical value of the ratio is approximately {final_value:.6f}.")

if __name__ == "__main__":
    solve_relativistic_aberration()