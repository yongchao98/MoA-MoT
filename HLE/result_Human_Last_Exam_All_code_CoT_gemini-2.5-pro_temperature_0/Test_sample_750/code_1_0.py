import numpy as np
import cmath
from scipy.optimize import fsolve

def solve_and_verify():
    """
    Finds and verifies solutions to the equation z*i = i^z.
    """
    print("Solving the equation z*i = i^z for complex numbers z.")
    print("-" * 50)

    # --- Part 1: Real Solutions ---
    # By analysis, we found two real solutions, z=1 and z=-1.
    # These solutions hold true for any integer branch k of the complex logarithm.
    
    print("Found a real solution: z = 1")
    z1 = 1
    # For z=1, the equation is (1)*i = i^(1)
    # We use k=0 for the principal value, but it holds for any k.
    log_i_k0 = 1j * (np.pi/2 + 2*0*np.pi)
    lhs1 = z1 * 1j
    rhs1 = cmath.exp(z1 * log_i_k0)
    print(f"Verifying the equation for z = {z1}:")
    print(f"({z1}) * i = {lhs1}")
    print(f"i^({z1}) = {rhs1}")
    print(f"The equation holds: {np.isclose(lhs1, rhs1)}")
    print("-" * 50)

    print("Found another real solution: z = -1")
    z2 = -1
    # For z=-1, the equation is (-1)*i = i^(-1)
    lhs2 = z2 * 1j
    rhs2 = cmath.exp(z2 * log_i_k0)
    print(f"Verifying the equation for z = {z2}:")
    print(f"({z2}) * i = {lhs2}")
    print(f"i^({z2}) = {rhs2}")
    print(f"The equation holds: {np.isclose(lhs2, rhs2)}")
    print("-" * 50)

    # --- Part 2: Purely Imaginary Solutions ---
    # A purely imaginary solution z = yi must satisfy -y = exp(-y * (pi/2 + 2*k*pi)).
    # Solutions for y exist only for k <= -1. We will find the solution for k = -1.
    
    k = -1
    print(f"Searching for a purely imaginary solution (using log branch k={k})...")
    
    # The equation to solve for y is: y + exp(y * (3*pi/2)) = 0
    # Let C = pi/2 + 2*k*pi = -3*pi/2
    # The equation is y + exp(-y*C) = 0
    C = np.pi/2 + 2*k*np.pi
    
    def imaginary_eq(y):
        return y[0] + np.exp(-y[0] * C)

    # We use fsolve to find the root. Analysis shows the root is negative.
    y_solution = fsolve(imaginary_eq, [-0.1])[0]
    z3 = y_solution * 1j
    
    print(f"Found a purely imaginary solution: z = {z3:.6f}")
    
    # For this z, the equation is (y*i)*i = i^(y*i)
    # We must use the corresponding log branch with k=-1.
    log_i_k_neg_1 = 1j * (np.pi/2 + 2*k*np.pi)
    lhs3 = z3 * 1j
    rhs3 = cmath.exp(z3 * log_i_k_neg_1)
    print(f"Verifying the equation for z = {z3:.6f}:")
    print(f"({z3:.6f}) * i = {lhs3:.6f}")
    print(f"i^({z3:.6f}) = {rhs3:.6f}")
    print(f"The equation holds: {np.isclose(lhs3, rhs3)}")
    print("-" * 50)
    
    print("Note: Infinitely many other complex solutions exist, which can be expressed")
    print("using the Lambert W function, but the ones above are the most fundamental.")

solve_and_verify()