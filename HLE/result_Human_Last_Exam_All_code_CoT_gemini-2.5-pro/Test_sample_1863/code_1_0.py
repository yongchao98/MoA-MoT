import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a
    quarter-sphere package, given certain constraints.
    """
    # 1. Define initial parameters
    package_diameter = 250.0

    # 2. Calculate the theoretical maximum diameter of the inscribed sphere.
    # The formula relating the package diameter (D) and the inscribed sphere
    # diameter (d) is: d = D / (1 + sqrt(3))
    sqrt_3 = math.sqrt(3)
    denominator = 1 + sqrt_3
    theoretical_diameter = package_diameter / denominator

    # 3. Apply the constraint that available sphere diameters are in increments of 0.01 cm.
    # We find the largest multiple of 0.01 that is less than or equal to the
    # theoretical maximum by truncating it to two decimal places.
    # math.floor(x * 100) / 100 achieves this.
    final_diameter = math.floor(theoretical_diameter * 100) / 100

    # 4. Print the explanation and the equation with its numerical values.
    print(f"The package is a quarter-sphere with a diameter of {package_diameter} cm.")
    print("The theoretical maximum diameter (d) of an inscribed sphere is found using the formula:")
    print("d = Package_Diameter / (1 + sqrt(3))\n")
    
    print("Plugging in the numbers:")
    print(f"d = {package_diameter} / (1 + {sqrt_3})")
    print(f"d = {package_diameter} / {denominator}")
    print(f"Theoretical d ≈ {theoretical_diameter} cm\n")

    print("Since available sphere diameters are in increments of 0.01 cm, we must use the largest possible value")
    print("that is not greater than the theoretical maximum. This means we truncate the value to two decimal places.\n")
    
    print(f"The maximum diameter of a sphere that can fit in the package is {final_diameter:.2f} cm.")

solve_sphere_in_quarter_sphere()
<<<91.50>>>