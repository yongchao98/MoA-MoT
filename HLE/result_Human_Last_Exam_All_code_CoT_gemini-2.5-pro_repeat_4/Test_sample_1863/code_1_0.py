import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Given diameter of the package's original sphere
    package_diameter = 250.0

    # 1. Calculate the radius of the quarter-sphere package
    package_radius = package_diameter / 2.0

    # 2. Formulate and solve the geometric relationship
    # The relationship between the package radius (R) and the inscribed sphere radius (r) is:
    # R = r * (sqrt(3) + 1)
    # So, the diameter d = 2 * r is:
    # d = 2 * R / (sqrt(3) + 1)
    # d = package_diameter / (sqrt(3) + 1)
    
    sqrt_3 = math.sqrt(3)
    theoretical_max_diameter = package_diameter / (sqrt_3 + 1)

    # 3. Apply the constraint of 0.01 cm increments
    # We need the largest diameter that is a multiple of 0.01 and <= theoretical_max_diameter.
    # This is equivalent to truncating the theoretical diameter to two decimal places.
    actual_max_diameter = math.floor(theoretical_max_diameter * 100) / 100

    # 4. Print the explanation and the result
    print("Step 1: Determine the radius of the quarter-sphere package.")
    print(f"The radius (R) is half the diameter: {package_diameter} cm / 2 = {package_radius} cm.\n")
    
    print("Step 2: Formulate the equation for the inscribed sphere's diameter (d).")
    print("The geometric relationship is R = r * (sqrt(3) + 1), where r is the inscribed sphere's radius.")
    print("This leads to the equation for the diameter d = 2 * r:\n")
    
    print("Final Equation:")
    print(f"d = (2 * R) / (sqrt(3) + 1)")
    print(f"d = (2 * {package_radius}) / ({sqrt_3} + 1)")
    print(f"d = {package_diameter} / {sqrt_3 + 1}")
    print(f"d \u2248 {theoretical_max_diameter} cm\n")

    print("Step 3: Apply the size constraint.")
    print("The available spheres have diameters in 0.01 cm increments.")
    print("Therefore, the maximum diameter is the theoretical value truncated to two decimal places.\n")

    print(f"The maximum diameter of a sphere that can fit in the package is: {actual_max_diameter} cm.")

solve_sphere_in_quarter_sphere()
<<<91.50>>>