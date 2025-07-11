import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # Diameter of the quarter-sphere package
    D_package = 250.0  # cm
    # Radius of the quarter-sphere package
    R_package = D_package / 2

    # The formula for the diameter 'd' of the inscribed sphere is d = D / (sqrt(3) + 1)
    # where D is the diameter of the large sphere.

    sqrt3 = math.sqrt(3)
    denominator = sqrt3 + 1
    
    # Calculate the theoretical maximum diameter
    d_theoretical = D_package / denominator

    # The available spheres have diameters in increments of 0.01 cm.
    # We must find the largest multiple of 0.01 that is less than or equal to the theoretical maximum.
    # This is achieved by flooring the value to two decimal places.
    d_final = math.floor(d_theoretical * 100) / 100

    print("Step 1: Define the relationship between the diameters.")
    print("The formula for the maximum theoretical diameter 'd' of a sphere inside a quarter-sphere of diameter 'D' is:")
    print("d = D / (sqrt(3) + 1)\n")
    
    print("Step 2: Substitute the known values into the equation.")
    print(f"d = {D_package} / ({sqrt3} + 1)")
    print(f"d = {D_package} / {denominator}")
    print(f"Theoretical maximum diameter d ≈ {d_theoretical} cm\n")

    print("Step 3: Apply the constraint of available sphere sizes (increments of 0.01 cm).")
    print("We must find the largest diameter that is a multiple of 0.01 and does not exceed the theoretical maximum.")
    print(f"Maximum allowed diameter = floor({d_theoretical:.6f} * 100) / 100")
    print(f"The maximum diameter of a sphere that can fit in the package is {d_final:.2f} cm.")

solve_sphere_in_quarter_sphere()
<<<91.50>>>