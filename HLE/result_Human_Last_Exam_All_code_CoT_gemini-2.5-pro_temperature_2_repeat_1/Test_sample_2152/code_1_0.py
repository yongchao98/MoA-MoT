import math

def calculate_ratio():
    """
    Calculates the ratio of differential cross-sections for scattering by a magnetic monopole
    and a magnetic dipole under the specified conditions.
    """
    
    # Scattering angle
    theta_rad = math.pi / 30
    
    # The final simplified ratio is R = 3 / (4 * sin^2(theta))
    # where theta is the scattering angle.
    
    # Calculate sin(theta)
    sin_theta = math.sin(theta_rad)
    
    # Calculate the ratio
    ratio = 3 / (4 * sin_theta**2)
    
    # Print the variables in the final equation
    print("The final expression for the ratio R is:")
    print("R = 3 / (4 * sin^2(theta))")
    print(f"where theta = pi/30 radians.")
    print("\nCalculation:")
    print(f"sin(pi/30) = {sin_theta}")
    print(f"sin^2(pi/30) = {sin_theta**2}")
    print(f"R = 3 / (4 * {sin_theta**2})")
    print(f"R = {ratio}")

if __name__ == '__main__':
    calculate_ratio()