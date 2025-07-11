import math

def solve():
    """
    Calculates the maximum diameter of a sphere that can fit inside a
    quarter-sphere package.
    """
    # 1. Define the package dimensions
    package_diameter = 250.0

    # 2. Derive the theoretical maximum diameter of the inscribed sphere.
    # The formula is d = D / (sqrt(3) + 1), where d is the inscribed sphere's
    # diameter and D is the package's original sphere diameter.
    sqrt_3 = math.sqrt(3)
    denominator = sqrt_3 + 1
    theoretical_diameter = package_diameter / denominator

    # 3. Adjust for the available sphere sizes.
    # The available diameters are in 0.01 cm increments, so we must round down
    # the theoretical maximum to find the largest sphere that can actually fit.
    max_diameter = math.floor(theoretical_diameter * 100) / 100.0

    # --- Output ---
    print("The formula for the maximum diameter (d) of a sphere that fits in a quarter-sphere of diameter (D) is:")
    print("d = D / (sqrt(3) + 1)\n")

    print("Plugging in the given values for the final equation:")
    # The prompt requires outputting each number in the final equation.
    print(f"d = {package_diameter} / (sqrt({3}) + {1})")
    print(f"d \u2248 {theoretical_diameter:.5f} cm\n")

    print("The available spheres have diameters in increments of 0.01 cm.")
    print("We must find the largest possible diameter, so we round the theoretical value down.\n")

    print(f"The maximum diameter of a sphere that can fit is {max_diameter:.2f} cm.")

solve()
<<<91.50>>>