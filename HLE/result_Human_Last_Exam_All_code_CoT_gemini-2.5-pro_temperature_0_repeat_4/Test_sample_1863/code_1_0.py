import math

def find_max_sphere_diameter():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # The diameter of the full sphere from which the package is made.
    package_diameter = 250.0

    # The relationship between the package diameter (D) and the maximum inner sphere diameter (d)
    # is given by the formula: d = D / (sqrt(3) + 1).

    # Calculate the theoretical maximum diameter.
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    theoretical_max_diameter = package_diameter / denominator

    # The available spheres have diameters in increments of 0.01 cm.
    # We must find the largest available diameter that is less than or equal to the theoretical maximum.
    # This is achieved by rounding down to two decimal places.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100

    # Print the calculation steps as requested.
    print("The equation for the theoretical maximum diameter (d) is: d = D / (sqrt(3) + 1)")
    print(f"Using the package diameter D = {package_diameter} cm, the calculation is:")
    print(f"d = {package_diameter} / ({sqrt_3} + 1) = {theoretical_max_diameter:.4f} cm")
    print("\nSince the sphere diameters are in increments of 0.01 cm, we round down.")
    print(f"The maximum diameter of a sphere that can fit in the package is: {final_diameter:.2f} cm")

find_max_sphere_diameter()