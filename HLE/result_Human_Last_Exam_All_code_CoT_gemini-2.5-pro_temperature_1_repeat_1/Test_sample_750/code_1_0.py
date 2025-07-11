import numpy as np
from scipy.special import lambertw

def solve_z_equation():
    """
    Finds and prints solutions to the equation z*i = i^z.
    """
    # 1. Real solutions
    print("Real solutions:")
    # For z = 1: 1*i = i and i^1 = i. This is a solution.
    z1 = 1
    print(f"z = {z1}")
    # For z = -1: -1*i = -i and i^-1 = 1/i = -i. This is a solution.
    z2 = -1
    print(f"z = {z2}")
    print("-" * 20)

    # 2. Purely imaginary solutions
    # These solutions are of the form z = iy, where x=0.
    # They exist for k = -m, where m is a positive integer (m >= 1).
    # The solutions for y can be expressed with the Lambert W function:
    # y_m = -W_0((4m-1)pi/2) / ((4m-1)pi/2)
    print("Purely imaginary solutions (first 5):")
    
    num_solutions_to_find = 5
    for m in range(1, num_solutions_to_find + 1):
        # Argument for the Lambert W function
        arg = (4 * m - 1) * np.pi / 2
        
        # The Lambert W function W_0(x) is real for x >= 0.
        # lambertw returns a complex number, so we take the real part.
        w_val = lambertw(arg, k=0).real
        
        y_m = -w_val / arg
        
        # The solution is z_m = i * y_m
        z_m = 1j * y_m
        
        # We need to output each number in the final equation,
        # which we interpret as printing the solution z.
        print(f"z = {z_m:.10f} (for k = {-m})")

    print("-" * 20)
    print("Note: An infinite number of other complex solutions (where x!=0 and y!=0) also exist.")

solve_z_equation()