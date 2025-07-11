import math

def find_max_sphere_diameter():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Given diameter of the quarter-sphere package in cm
    package_diameter = 250.0

    # The relationship between the package diameter (D) and the max inner sphere diameter (d) is:
    # d = D / (sqrt(3) + 1)

    # --- Calculation Steps ---
    
    # 1. Calculate the values for the equation
    sqrt3_val = math.sqrt(3)
    denominator = sqrt3_val + 1
    
    # 2. Calculate the theoretical maximum diameter
    max_diameter_theoretical = package_diameter / denominator
    
    # 3. Adjust for the 0.01 cm increment constraint by flooring to two decimal places
    final_diameter = math.floor(max_diameter_theoretical * 100) / 100

    # --- Output the result ---

    print("The formula for the maximum inner sphere diameter (d) is: d = D / (sqrt(3) + 1)")
    print(f"Where the package diameter D = {package_diameter} cm.")
    print("\nThe final equation with the numbers plugged in is:")
    # Using f-string to display each number in the equation
    print(f"{package_diameter} / ( {sqrt3_val} + 1 ) = {max_diameter_theoretical}")
    
    print("\nSince the spheres are available in increments of 0.01 cm, we take the largest possible size that fits.")
    print(f"Maximum allowed diameter = {final_diameter} cm")

if __name__ == '__main__':
    find_max_sphere_diameter()