import math

def calculate_cross_section_ratio():
    """
    Calculates the ratio of differential cross-sections for scattering
    off a magnetic monopole versus a magnetic dipole.
    """
    # Given parameters
    # Ratio of speed of light to particle speed, c/v
    c_v_ratio = 100
    # Scattering angle in radians
    theta_rad = math.pi / 30

    # The derived formula for the ratio R is:
    # R = (3 / 40000) * (c/v)^2 / sin(theta)^2
    # This simplifies to R = (3 * 100^2) / (40000 * sin(theta)^2)
    # R = 30000 / (40000 * sin(theta)^2) = 3 / (4 * sin(theta)^2)

    # The numbers in the final simplified equation are 3, 4, pi, and 30.
    numerator = 3
    denominator_const = 4
    angle_denominator = 30

    # Calculate the value of sin(theta)^2
    sin_theta = math.sin(math.pi / angle_denominator)
    sin_theta_sq = sin_theta**2

    # Calculate the final ratio
    ratio = numerator / (denominator_const * sin_theta_sq)

    # Print the explanation and the result
    print("The final equation for the ratio of differential cross-sections (monopole/dipole) is:")
    print(f"Ratio = {numerator} / ({denominator_const} * sin(pi/{angle_denominator})^2)")
    print("\nBreaking down the calculation:")
    print(f"The value of sin(pi/{angle_denominator}) is {sin_theta}")
    print(f"The value of sin(pi/{angle_denominator})^2 is {sin_theta_sq}")
    print(f"The equation becomes: Ratio = {numerator} / ({denominator_const} * {sin_theta_sq})")
    print(f"Final calculated ratio: {ratio}")

calculate_cross_section_ratio()