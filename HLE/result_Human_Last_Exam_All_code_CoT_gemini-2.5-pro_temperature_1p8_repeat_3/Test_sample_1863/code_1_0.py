import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # 1. Define the properties of the package.
    package_diameter = 250.0  # cm

    # 2. Formulate and solve the geometric relationship.
    # The relationship between the package diameter (D) and the maximum
    # inner sphere diameter (d) is d = D / (sqrt(3) + 1).
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    theoretical_diameter = package_diameter / denominator

    # 3. Apply the manufacturing constraint.
    # The sphere's diameter must be a multiple of 0.01 cm.
    # We find the largest such diameter that is <= theoretical_diameter
    # by flooring the value to 2 decimal places.
    final_diameter = math.floor(theoretical_diameter * 100) / 100

    # 4. Print the explanation and result.
    print("The maximum diameter (d) of a sphere that fits in a quarter-sphere package (diameter D) is found with the formula:")
    print("d = D / (sqrt(3) + 1)\n")
    print("Substituting the given values into the equation:")
    print(f"d = {package_diameter} / (sqrt({3}) + 1)")
    print(f"d = {package_diameter} / ({sqrt_3} + 1)")
    print(f"d = {package_diameter} / {denominator}")
    print(f"d ≈ {theoretical_diameter} cm\n")

    print("The available spheres have diameters in increments of 0.01 cm.")
    print("Therefore, the maximum diameter must be rounded down to the nearest 0.01 cm.\n")

    print(f"Final Answer: The maximum diameter of the sphere is {final_diameter} cm.")


solve_sphere_in_quarter_sphere()
<<<91.50>>>