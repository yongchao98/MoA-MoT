import math

def calculate_jet_speed(diameter_mm):
    """
    Calculates the speed of a gas jet from a bursting bubble using a model
    based on Laplace pressure and Bernoulli's principle.

    Args:
        diameter_mm (float): The diameter of the bubble in millimeters.

    Returns:
        float: The calculated speed of the gas jet in m/s.
    """
    # Physical constants
    sigma = 0.072  # Surface tension of water (N/m)
    rho_air = 1.225  # Density of air (kg/m^3)

    # Convert diameter to radius in meters
    diameter_m = diameter_mm / 1000.0
    radius_m = diameter_m / 2.0

    # Calculate jet speed using the formula v = sqrt(4*sigma / (R * rho_air))
    jet_speed = math.sqrt((4 * sigma) / (radius_m * rho_air))

    print(f"For a bubble diameter of {diameter_mm} mm (radius = {radius_m:.4f} m):")
    print(f"v = sqrt( (4 * {sigma}) / ({radius_m:.4f} * {rho_air}) )")
    print(f"v = {jet_speed:.2f} m/s")
    print("-" * 20)
    return jet_speed

# --- Main execution ---
# Calculate for a 2 mm diameter bubble
speed1 = calculate_jet_speed(2)

# Calculate for a 2 cm (20 mm) diameter bubble
speed2 = calculate_jet_speed(20)

print(f"The calculated speeds are approximately {round(speed1)} m/s and {round(speed2, 1)} m/s.")
print("Comparing these results with the answer choices, the first value (15 m/s) is matched.")
