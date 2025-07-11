import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a
    quarter-sphere package.
    """
    # Diameter of the quarter-sphere package in cm.
    package_diameter = 250.0

    # The formula for the diameter of the inscribed sphere (d) in terms of the
    # package's diameter (D) is: d = D / (sqrt(2) + 1)
    
    # Calculate the components of the formula
    sqrt_2 = math.sqrt(2)
    denominator = sqrt_2 + 1
    
    # Calculate the final diameter of the inscribed sphere
    sphere_diameter = package_diameter / denominator
    
    print("The formula to find the maximum diameter of the inscribed sphere (d) is:")
    print("d = D / (sqrt(2) + 1)")
    print(f"Where D (package diameter) = {package_diameter} cm.\n")
    print("Substituting the values into the equation:")
    print(f"d = {package_diameter} / ({sqrt_2} + 1)")
    print(f"d = {sphere_diameter}\n")
    print("Therefore, the maximum diameter of a sphere that can fit in the package is:")
    print(f"{sphere_diameter} cm")

solve_sphere_in_quarter_sphere()