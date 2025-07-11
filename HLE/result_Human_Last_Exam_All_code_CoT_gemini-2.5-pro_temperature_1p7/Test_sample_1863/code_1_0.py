import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Diameter of the quarter-sphere package in cm
    package_diameter = 250.0

    # The formula for the theoretical maximum diameter (d) of the inscribed sphere is
    # d = D / (sqrt(3) + 1), where D is the diameter of the package sphere.

    # Calculate the necessary components of the formula
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1

    # Calculate the theoretical maximum diameter
    theoretical_max_diameter = package_diameter / denominator

    # The available spheres have diameters in 0.01 cm increments.
    # We must find the largest available sphere that can fit, so we take the floor
    # of the theoretical maximum to two decimal places.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100

    # Output the explanation and the final equation with its numbers
    print("The formula to find the maximum theoretical diameter (d) is: d = D / (sqrt(3) + 1)")
    print("\nPlugging in the package diameter D = 250.0 cm:")
    print(f"d = {package_diameter} / ({sqrt_3} + 1)")
    print(f"d = {package_diameter} / {denominator}")
    print(f"Theoretical maximum diameter = {theoretical_max_diameter} cm")

    print("\nSince the available spheres are in increments of 0.01 cm, the actual maximum diameter is the theoretical value floored to 2 decimal places.")
    print(f"Maximum diameter = {final_diameter:.2f} cm")

solve_sphere_in_quarter_sphere()