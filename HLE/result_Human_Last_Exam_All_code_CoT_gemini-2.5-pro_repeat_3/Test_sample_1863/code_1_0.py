import math

def find_max_sphere_diameter():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # 1. Define the parameters from the problem description
    package_diameter = 250.0
    diameter_increment = 0.01

    # 2. Derive the formula for the theoretical maximum diameter.
    # The formula for the diameter 'd' of the largest sphere that can be inscribed
    # in a quarter-sphere of diameter 'D' is d = D * (sqrt(2) - 1).
    # This is derived from the geometric constraint that the inscribed sphere must be
    # tangent to the two flat faces and the curved face of the quarter-sphere.
    
    # 3. Calculate the theoretical maximum diameter
    sqrt_2 = math.sqrt(2)
    theoretical_max_diameter = package_diameter * (sqrt_2 - 1)

    # 4. Find the largest available diameter given the 0.01 cm increment.
    # The sphere must fit entirely, so its diameter must be <= the theoretical max.
    # We find the largest multiple of the increment that meets this condition.
    final_diameter = math.floor(theoretical_max_diameter / diameter_increment) * diameter_increment
    
    # 5. Output the explanation and results, showing the equation and numbers used.
    print("The formula for the theoretical maximum diameter (d) of a sphere inscribed in a quarter-sphere of diameter D is:")
    print("d = D * (sqrt(2) - 1)")
    print("\nPlugging in the given numbers:")
    print(f"d = {package_diameter} * (math.sqrt(2) - 1)")
    print(f"d = {package_diameter} * ({sqrt_2} - 1)")
    print(f"Theoretical maximum diameter = {theoretical_max_diameter} cm")

    print("\nSince sphere diameters are available in increments of 0.01 cm, we must find the largest multiple of 0.01 that is not greater than the theoretical maximum.")
    print(f"Maximum practical diameter = floor({theoretical_max_diameter} / {diameter_increment}) * {diameter_increment}")
    print(f"Maximum practical diameter = {final_diameter} cm")

if __name__ == '__main__':
    find_max_sphere_diameter()