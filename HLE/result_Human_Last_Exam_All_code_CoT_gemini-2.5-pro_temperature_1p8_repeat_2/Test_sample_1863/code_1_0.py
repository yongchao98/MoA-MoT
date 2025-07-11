import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit in a quarter-sphere package.
    """
    # 1. Define problem parameters
    package_diameter = 250.0
    
    # 2. Explain the formula
    print("To find the maximum diameter (d) of a sphere that fits in a quarter-sphere package of diameter (D),")
    print("we use the geometric relationship derived from placing the sphere tangent to all four surfaces.")
    print("This gives the formula:")
    print("d = D / (1 + sqrt(3))")
    print("-" * 30)

    # 3. Calculate the theoretical maximum diameter
    print(f"Given the package diameter D = {package_diameter} cm.")
    print("\nSubstitute the values into the formula:")
    
    # Show the equation with numbers
    sqrt_3 = math.sqrt(3)
    denominator = 1 + sqrt_3
    max_theoretical_diameter = package_diameter / denominator
    
    print(f"d = {package_diameter} / (1 + √3)")
    print(f"d = {package_diameter} / (1 + {sqrt_3:.5f})")
    print(f"d = {package_diameter} / {denominator:.5f}")
    print(f"d ≈ {max_theoretical_diameter:.5f} cm")
    print("-" * 30)

    # 4. Adjust for the 0.01 cm increment constraint
    print("The available spheres have diameters in increments of 0.01 cm.")
    print("We must find the largest multiple of 0.01 cm that does not exceed the theoretical maximum.")
    
    # Truncate to two decimal places
    final_diameter = math.floor(max_theoretical_diameter * 100) / 100
    
    print(f"\nFinal Diameter = floor({max_theoretical_diameter:.5f} * 100) / 100")
    print(f"Final Diameter = {final_diameter:.2f} cm")

solve_sphere_in_quarter_sphere()