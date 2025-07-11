import math

def find_max_sphere_diameter():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # The diameter of the full sphere from which the quarter-sphere package is made.
    package_diameter = 250.0

    # The formula for the diameter (d) of the largest sphere that can be inscribed
    # in a quarter-sphere of original diameter (D) is derived from geometric constraints:
    # d = D / (sqrt(3) + 1)

    # Calculate the values needed for the equation.
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    
    # Calculate the theoretical maximum diameter.
    theoretical_max_diameter = package_diameter / denominator

    # The available spheres have diameters in increments of 0.01 cm.
    # We must find the largest multiple of 0.01 that is less than or equal to the theoretical maximum.
    # This is achieved by truncating the value to two decimal places.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100

    # Print the final equation with each number, as requested.
    print("The equation to find the maximum diameter (d) is based on the package's original diameter (D):")
    print("d = D / (sqrt(3) + 1)")
    print("\nPlugging in the numbers:")
    print(f"d = {package_diameter} / ({sqrt_3} + 1)")
    print(f"d = {package_diameter} / {denominator}")
    print(f"d \u2248 {theoretical_max_diameter:.6f} cm")
    
    print("\nSince the available sphere diameters are in increments of 0.01 cm, we must use the largest available size that fits.")
    print(f"The maximum diameter of such a sphere is: {final_diameter:.2f} cm")

find_max_sphere_diameter()