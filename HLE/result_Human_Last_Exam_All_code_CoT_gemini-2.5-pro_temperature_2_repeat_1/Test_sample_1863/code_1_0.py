import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Given parameters
    package_diameter = 250.0
    
    # --- Step 1: Calculate the radius of the package ---
    package_radius = package_diameter / 2
    
    # --- Step 2: Explain the geometric derivation ---
    print(f"The package is a quarter-sphere with a diameter of {package_diameter} cm.")
    print(f"Its radius (R) is therefore {package_diameter} / 2 = {package_radius} cm.")
    print("\nTo find the maximum diameter (d) of a sphere that can fit inside, we use the geometric relationship:")
    print("d = (Package Diameter) / (1 + sqrt(3))")
    print("\nLet's calculate this value:")
    
    # --- Step 3: Calculate the maximum theoretical diameter ---
    sqrt_3 = math.sqrt(3)
    denominator = 1 + sqrt_3
    max_theoretical_diameter = package_diameter / denominator
    
    # --- Step 4: Show the equation with numbers plugged in ---
    print(f"d = {package_diameter} / (1 + {sqrt_3:.5f})")
    print(f"d = {package_diameter} / {denominator:.5f}")
    print(f"d \u2248 {max_theoretical_diameter:.5f} cm")
    
    # --- Step 5: Apply the 0.01 cm increment constraint ---
    # We round down to the nearest 0.01 cm because the sphere must fit entirely inside.
    max_practical_diameter = math.floor(max_theoretical_diameter * 100) / 100
    
    # --- Step 6: Print the final answer ---
    print("\nSince the available spheres have diameters in increments of 0.01 cm, we must round down to find the largest one that fits.")
    print(f"The maximum diameter of a sphere that can be packed is: {max_practical_diameter:.2f} cm")

solve_sphere_in_quarter_sphere()