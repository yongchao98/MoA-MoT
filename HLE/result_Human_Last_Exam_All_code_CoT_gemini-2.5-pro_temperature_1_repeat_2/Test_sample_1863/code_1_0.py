import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a
    quarter-sphere package, considering diameter increments of 0.01 cm.
    """
    # Diameter of the quarter-sphere package in cm.
    package_diameter = 250

    # The formula for the maximum theoretical diameter (d) of the inscribed sphere is:
    # d = D / (sqrt(3) + 1), where D is the diameter of the package.
    
    # Calculate the components of the equation
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    
    # Calculate the theoretical maximum diameter
    max_theoretical_diameter = package_diameter / denominator
    
    # The available spheres have diameters in increments of 0.01 cm.
    # We must round the theoretical diameter down to the nearest 0.01 cm.
    final_diameter = math.floor(max_theoretical_diameter * 100) / 100

    # Output the steps of the final equation as requested
    print("The formula for the maximum diameter (d) is: d = D / (sqrt(3) + 1)")
    print(f"Using D = {package_diameter}, the equation is:")
    print(f"d = {package_diameter} / (math.sqrt(3) + 1)")
    print(f"d = {package_diameter} / ({sqrt_3} + 1)")
    print(f"d = {package_diameter} / {denominator}")
    print(f"Theoretical maximum d ≈ {max_theoretical_diameter} cm")
    
    print("\nRounding down to the nearest available 0.01 cm increment:")
    print(f"The maximum diameter of a sphere that can fit is {final_diameter:.2f} cm.")

solve_sphere_in_quarter_sphere()