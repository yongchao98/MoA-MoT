import math

def calculate_angular_distance():
    """
    Calculates the angular distance between the two stars.

    The problem provides several pieces of information about two stars, A and B,
    whose apparent positions change only due to Earth's axial precession.
    - Precession period T = 26000 years.
    - Axial tilt (obliquity) ε = 23.5 degrees.
    - Star A was last on the celestial equator 3000 years ago.
    - Star B will first be on the celestial equator 10000 years from now.
    - The time between these two equator crossings is 10000 - (-3000) = 13000 years.
    - This interval is exactly half the precession period (26000 / 2 = 13000 years).
    - During this 13000-year interval, the stars are on the same side of the celestial equator.
    - A future swap of equatorial coordinates implies a deep symmetry between the stars.

    This setup strongly implies a specific symmetrical configuration. The simplest configuration
    that satisfies all conditions is one where the two stars are located on the path traced by
    the celestial pole, but separated by 180 degrees of precession.
    This corresponds to the stars being at opposite points on the circle of radius ε around the
    ecliptic pole. The angular distance between two such points is the diameter of this circle,
    which is twice the radius.

    Therefore, the angular distance `d` between the two stars is 2 * ε.
    """
    
    # Earth's axial tilt in degrees
    axial_tilt_epsilon = 23.5
    
    # The angular distance is twice the axial tilt
    angular_distance = 2 * axial_tilt_epsilon
    
    # We print the final equation
    print(f"The angular distance is calculated as twice the Earth's axial tilt.")
    print(f"Distance = 2 * {axial_tilt_epsilon} degrees")
    print(f"Distance = {angular_distance} degrees")

calculate_angular_distance()