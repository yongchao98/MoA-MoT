import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that fits in a quarter-sphere package.
    """
    # 1. Define the geometry of the package
    package_diameter = 250.0
    package_radius = package_diameter / 2

    # 2. Use the derived formula to find the theoretical maximum diameter of the inner sphere.
    # The formula is: d = 2 * R / (sqrt(3) + 1), where R is the package radius.
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    theoretical_diameter = (2 * package_radius) / denominator

    # 3. Adjust for the available sphere sizes (increments of 0.01 cm).
    # We need the largest diameter <= theoretical_diameter. This is achieved by
    # truncating (or flooring) the result to two decimal places.
    final_diameter = math.floor(theoretical_diameter * 100) / 100

    # 4. Print the steps and the final answer.
    print("Step 1: Determine the radius of the quarter-sphere package.")
    print(f"Package Diameter = {package_diameter} cm")
    print(f"Package Radius (R) = {package_diameter} / 2 = {package_radius} cm")
    print("-" * 30)

    print("Step 2: Use the formula for the maximum diameter (d) of an inscribed sphere.")
    print("Formula: d = (2 * R) / (sqrt(3) + 1)")
    print("Substituting the values:")
    print(f"d = (2 * {package_radius}) / ({sqrt_3} + 1)")
    print(f"d = {2 * package_radius} / {denominator}")
    print(f"Theoretical maximum diameter = {theoretical_diameter} cm")
    print("-" * 30)

    print("Step 3: Find the largest available diameter.")
    print("The available spheres have diameters in increments of 0.01 cm.")
    print("We must round down the theoretical diameter to the nearest 0.01 cm.")
    print(f"Final Answer: The maximum diameter of a sphere that can fit is {final_diameter:.2f} cm.")

if __name__ == "__main__":
    solve_sphere_in_quarter_sphere()