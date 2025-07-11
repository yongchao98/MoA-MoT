import math

def calculate_angular_distance():
    """
    Calculates the angular distance between two stars based on precession data.

    The problem's constraints lead to a unique scenario where the stars'
    ecliptic latitude (β) must be equal to the Earth's axial tilt (ε).
    The angular distance (d) between two stars at the same ecliptic latitude β
    but opposite ecliptic longitudes is given by the formula d = 180 - 2β.
    """
    
    # Earth's axial tilt in degrees, given in the problem
    axial_tilt_epsilon = 23.5
    
    # The ecliptic latitude of the stars must be equal to the axial tilt
    beta = axial_tilt_epsilon
    
    # Calculate the angular distance in degrees
    # The formula is derived from the spherical law of cosines for two points
    # with coordinates (lambda, beta) and (lambda + 180, beta).
    # cos(d) = sin(beta)*sin(beta) + cos(beta)*cos(beta)*cos(180)
    # cos(d) = sin^2(beta) - cos^2(beta) = -cos(2*beta)
    # d = arccos(-cos(2*beta)) = 180 - 2*beta
    angular_distance = 180.0 - 2 * beta
    
    # Print the equation with the final numbers
    print(f"The angular distance 'd' is calculated using the formula: d = 180 - 2 * \u03B2")
    print(f"Given the constraints, the ecliptic latitude \u03B2 must equal the Earth's axial tilt \u03B5.")
    print(f"d = 180 - 2 * {beta}")
    
    # Print the final result
    print(f"The angular distance between the two stars is: {angular_distance} degrees")

calculate_angular_distance()