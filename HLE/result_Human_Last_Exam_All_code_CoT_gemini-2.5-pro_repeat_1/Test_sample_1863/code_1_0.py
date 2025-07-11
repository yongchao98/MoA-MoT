import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a
    quarter-sphere package.
    """
    # 1. Define the parameters of the package
    package_diameter = 250.0
    package_radius = package_diameter / 2.0

    # 2. Derive the formula for the inner sphere's diameter
    # Let R be the radius of the large quarter-sphere (package).
    # Let r be the radius of the small sphere to be fitted inside.
    # The center of the small sphere will be at (r, r, r).
    # The distance from the origin to this center is sqrt(3)*r.
    # The tangency condition is: sqrt(3)*r + r = R
    # Solving for the diameter d = 2r:
    # r * (sqrt(3) + 1) = R
    # d/2 * (sqrt(3) + 1) = R
    # d = 2 * R / (sqrt(3) + 1)
    # A simplified version of the formula is d = R * (sqrt(3) - 1)

    # 3. Calculate the theoretical maximum diameter
    sqrt_3 = math.sqrt(3)
    theoretical_max_diameter = package_radius * (sqrt_3 - 1)

    # 4. Adjust for the 0.01 cm increment constraint
    # We must find the largest multiple of 0.01 that is less than or
    # equal to the theoretical maximum. This is done by truncating
    # the value to two decimal places.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100

    # 5. Print the results showing the steps
    print(f"The radius of the quarter-sphere package (R) is {package_diameter} / 2 = {package_radius} cm.")
    print("The formula for the maximum theoretical diameter (d) of the inscribed sphere is: d = R * (sqrt(3) - 1)")
    print("Plugging in the numbers:")
    print(f"d = {package_radius} * ({sqrt_3} - 1)")
    print(f"d = {package_radius} * {sqrt_3 - 1}")
    print(f"Theoretical maximum diameter ≈ {theoretical_max_diameter} cm.")
    print("\nSince available diameters are in 0.01 cm steps, we take the largest possible value by truncating to two decimal places.")
    print(f"The maximum diameter of a sphere that can fit in the package is {final_diameter} cm.")

solve_sphere_in_quarter_sphere()