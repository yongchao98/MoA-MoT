import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # The diameter of the original sphere from which the quarter-sphere is made.
    package_diameter = 250.0  # in cm

    # The radius of the quarter-sphere package.
    R = package_diameter / 2.0

    # The relationship between the radius of the quarter-sphere (R) and the radius
    # of the maximum inscribed sphere (r) is: r * (sqrt(2) + 1) = R.
    # Therefore, the diameter of the inscribed sphere (d = 2r) is related to the
    # diameter of the package (D = 2R) by: d = D / (sqrt(2) + 1).

    # Calculate the values for the equation.
    sqrt_2 = math.sqrt(2)
    denominator = sqrt_2 + 1
    max_diameter = package_diameter / denominator

    # Print the explanation and the final equation with all numbers.
    print(f"The diameter of the quarter-sphere package is D = {package_diameter} cm.")
    print("The maximum diameter 'd' of a sphere that can fit inside is calculated by the formula: d = D / (sqrt(2) + 1)")
    print("\nFinal Equation:")
    print(f"d = {package_diameter} / (sqrt(2) + 1)")
    print(f"d = {package_diameter} / ({sqrt_2} + 1)")
    print(f"d = {package_diameter} / {denominator}")
    print(f"d = {max_diameter} cm")

solve_sphere_in_quarter_sphere()

# The problem asks for the maximum diameter of such a sphere.
# The calculation above gives the result.
# The final answer is the value of max_diameter.
# max_diameter = 250 / (math.sqrt(2) + 1) which is approximately 103.55339
print(f"\n<<<{250 / (math.sqrt(2) + 1)}>>>")