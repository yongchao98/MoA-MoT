import numpy as np
from scipy.special import lambertw
import cmath

def solve_and_print_solutions():
    """
    Solves the equation z*i = i^z and prints the solutions.
    """
    print("The equation z * i = i^z has two real solutions and an infinite set of purely imaginary solutions.")
    print("-" * 40)

    # --- Real Solutions ---
    print("Real solutions:")
    # Solution z=1
    z1 = 1
    lhs1_val = z1 * 1j
    rhs1_val = 1j**z1
    print(f"Solution 1: z = {z1}")
    print(f"Checking the equation: ({z1}) * i = i^({z1})")
    print(f"Result: {lhs1_val} = {rhs1_val}")
    print()

    # Solution z=-1
    z2 = -1
    lhs2_val = z2 * 1j
    rhs2_val = 1j**z2
    print(f"Solution 2: z = {z2}")
    print(f"Checking the equation: ({z2}) * i = i^({z2})")
    print(f"Result: {lhs2_val:.6f} = {rhs2_val:.6f}")
    print("-" * 40)

    # --- Imaginary Solutions ---
    print("Purely imaginary solutions:")
    print("These solutions z = i*y exist for each integer k <= -1.")
    print("y_k is found by solving -y = exp(-y * C_k), where C_k = pi/2 + 2*k*pi.")
    print("Using the Lambert W function, y_k = W(-C_k) / C_k.")
    print("The first 5 imaginary solutions (for k = -1 to -5) are:")

    for k in range(-1, -6, -1):
        # C_k = pi/2 + 2k*pi
        C_k = np.pi / 2 + 2 * k * np.pi
        
        # Argument for the Lambert W function is -C_k. For k <= -1, -C_k is positive.
        arg_w = -C_k
        
        # For a positive real argument, the principal branch W_0 is real.
        # The lambertw function returns a complex number, so we take the real part.
        w_val = lambertw(arg_w).real
        
        # Calculate y_k
        y_k = w_val / C_k
        
        # The imaginary solution z_k
        z_k = y_k * 1j
        
        # For imaginary solutions, the original equation z*i = i^z simplifies.
        # LHS = (i*y)*i = -y
        # RHS = i^(i*y) = exp(i*y * log(i)) = exp(i*y * i*C_k) = exp(-y*C_k)
        # So, the equation we solve is -y_k = exp(-y_k * C_k)
        final_eq_lhs = -y_k
        final_eq_rhs = np.exp(-y_k * C_k)

        print(f"\nSolution for k = {k}:")
        print(f"z = {z_k:.8f}")
        print(f"The equation becomes -y = exp(-y*C_k), with y = {y_k:.8f}")
        print(f"Checking the simplified equation: {final_eq_lhs:.8f} = {final_eq_rhs:.8f}")

# Execute the function to find and print solutions
solve_and_print_solutions()