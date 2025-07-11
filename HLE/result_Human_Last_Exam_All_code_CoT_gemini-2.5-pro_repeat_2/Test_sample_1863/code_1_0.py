import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that fits inside a quarter-sphere package.
    """
    # Step 1: Define the dimensions of the package
    package_diameter = 250.0
    package_radius = package_diameter / 2.0

    # Step 2: Explain the geometric relationship
    # The diameter 'd' of the inner sphere is related to the radius 'R' of the
    # quarter-sphere by the formula: d = (2 * R) / (sqrt(3) + 1)
    # This simplifies to: d = package_diameter / (sqrt(3) + 1)
    print("The goal is to find the maximum diameter of a sphere that can fit inside a quarter-sphere package.")
    print(f"The package has a diameter of {package_diameter} cm, so its radius (R) is {package_radius} cm.")
    print("\nThe formula for the maximum diameter (d) of the inscribed sphere is:")
    print("d = (2 * R) / (sqrt(3) + 1)")
    
    # Step 3: Substitute the values and show the equation
    # We use the package diameter directly for a simpler equation: d = package_diameter / (sqrt(3) + 1)
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    
    print("\nSubstituting the package diameter:")
    print(f"d = {package_diameter} / (sqrt(3) + 1)")
    print(f"d = {package_diameter} / ({sqrt_3:.4f} + 1)")
    print(f"d = {package_diameter} / {denominator:.4f}")

    # Step 4: Calculate the theoretical maximum diameter
    theoretical_max_diameter = package_diameter / denominator
    print(f"d \u2248 {theoretical_max_diameter:.4f} cm")

    # Step 5: Adjust for the available sphere sizes (increments of 0.01 cm)
    # We must round down to the nearest 0.01 cm to ensure it fits.
    # This is equivalent to truncating to two decimal places.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100
    
    print("\nSince available spheres have diameters in 0.01 cm increments, we round down to the nearest valid size.")
    print(f"The maximum diameter of a sphere that can fit in the package is {final_diameter:.2f} cm.")

solve_sphere_in_quarter_sphere()