import math

def calculate_star_distance():
    """
    Calculates the angular distance between two stars based on the given orbital mechanics problem.
    """
    # The Earth's axial tilt (obliquity of the ecliptic), epsilon, in degrees.
    axial_tilt_deg = 23.5

    # Based on the problem's symmetry and timing conditions, it is derived that the
    # magnitude of the ecliptic latitude for both stars, |beta|, must be equal
    # to the Earth's axial tilt, epsilon.
    beta_deg = axial_tilt_deg

    # The angular distance 'd' between two stars with the same ecliptic latitude 'beta'
    # and opposite ecliptic longitudes (separated by 180 degrees) is given by the formula:
    # d = 180 - 2 * |beta|
    constant_angle = 180
    distance_deg = constant_angle - 2 * beta_deg

    # Output the final calculation, showing each number in the equation.
    print("The angular distance, d, is calculated using the formula derived from the stars' relative positions:")
    print(f"d = {constant_angle} - 2 * |beta|")
    print(f"With |beta| being equal to the axial tilt, {beta_deg} degrees, the calculation is:")
    print(f"d = {constant_angle} - 2 * {beta_deg} = {distance_deg}")
    print(f"\nThe final angular distance between the two stars is {distance_deg} degrees.")

calculate_star_distance()