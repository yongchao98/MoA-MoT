import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Define the parameters from the problem description
    package_diameter = 250.0
    diameter_increment = 0.01

    # Step 1: Establish the geometric relationship to find the max theoretical diameter.
    # Let D be the diameter of the full sphere from which the quarter-sphere is cut.
    # Let d be the diameter of the smaller sphere to be fitted inside.
    # For the small sphere to be maximum, it must be tangent to the three flat faces
    # and the one curved face of the quarter-sphere.
    # This geometric constraint leads to the following equation:
    # d * (sqrt(3) + 1) = D

    # Step 2: Solve the equation for the theoretical maximum diameter 'd_max'.
    # We substitute the given package diameter D = 250 cm.
    print("The formula for the maximum theoretical diameter 'd_max' is derived from the geometry:")
    print("d_max * (sqrt(3) + 1) = D")
    print("\nSolving for d_max:")
    print("d_max = D / (sqrt(3) + 1)")
    
    sqrt_3 = math.sqrt(3)
    theoretical_max_diameter = package_diameter / (sqrt_3 + 1)

    print("\nSubstituting the numbers into the equation:")
    # As requested, outputting each number in the final equation.
    print(f"d_max = {package_diameter} / ( {sqrt_3:.5f} + 1.0 )")
    print(f"d_max = {package_diameter} / {sqrt_3 + 1.0:.5f}")
    print(f"The calculated maximum theoretical diameter is: {theoretical_max_diameter:.5f} cm.")

    # Step 3: Account for the available sphere sizes.
    # The diameters are available in increments of 0.01 cm. We need the largest
    # multiple of 0.01 that is less than or equal to the theoretical maximum.
    # This is equivalent to rounding down to two decimal places.
    
    num_increments = math.floor(theoretical_max_diameter / diameter_increment)
    final_diameter = num_increments * diameter_increment

    print(f"\nSince sphere diameters are available in increments of {diameter_increment} cm, we find the largest valid diameter:")
    print(f"The final maximum diameter is floor({theoretical_max_diameter:.5f} / {diameter_increment}) * {diameter_increment}")
    print(f"\nThe maximum diameter of a sphere that can fit in the package is: {final_diameter:.2f} cm.")

solve_sphere_in_quarter_sphere()