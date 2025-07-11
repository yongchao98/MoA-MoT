import numpy as np
from scipy.special import lambertw
import cmath

def solve_and_verify_complex_equation():
    """
    Solves for z in the complex equation z*i = i^z, prints, and verifies the solutions.
    """
    print("The equation to solve is z * i = i^z for z in complex numbers.")
    print("This equation has two real solutions and an infinite family of purely imaginary solutions.\n")

    # --- Real Solutions ---
    print("--- Real Solutions ---")
    
    # Solution 1: z = 1
    z1 = 1 + 0j
    lhs1 = z1 * 1j
    # For z=1, i^1 = exp(1 * log(i)) = exp(i*(pi/2 + 2k*pi)) = i, for any integer k.
    # cmath.pow uses the principal value for the logarithm.
    rhs1 = cmath.pow(1j, z1) 
    print(f"Solution 1: z = {z1.real}")
    print(f"The equation is: ({z1.real}) * i = i^({z1.real})")
    print(f"Verification: LHS = {lhs1:.6f}, RHS = {rhs1:.6f}. Match: {np.isclose(lhs1, rhs1)}\n")

    # Solution 2: z = -1
    z2 = -1 + 0j
    lhs2 = z2 * 1j
    # For z=-1, i^-1 = exp(-1 * log(i)) = exp(-i*(pi/2 + 2k*pi)) = -i, for any integer k.
    rhs2 = cmath.pow(1j, z2)
    print(f"Solution 2: z = {z2.real}")
    print(f"The equation is: ({z2.real}) * i = i^({z2.real})")
    print(f"Verification: LHS = {lhs2:.6f}, RHS = {rhs2:.6f}. Match: {np.isclose(lhs2, rhs2)}\n")

    # --- Purely Imaginary Solutions ---
    print("--- Purely Imaginary Solutions ---")
    print("An infinite family of purely imaginary solutions exists, one for each negative integer k (k <= -1).")
    print("These solutions arise from the different branches of the complex logarithm used to define i^z.")
    print("i^z = exp(z * (log|i| + i*arg(i))) = exp(z * i * (pi/2 + 2*k*pi))\n")
    print("Calculating the first 5 solutions for k = -1 down to -5:\n")

    for k in range(-1, -6, -1):
        # For a given k <= -1, we solve -y = exp(-y * (pi/2 + 2*k*pi)) for z=iy.
        # Let alpha = -(pi/2 + 2*k*pi) = pi*(-2*k - 1/2). Since k <= -1, alpha > 0.
        alpha = np.pi * (-2 * k - 0.5)
        
        # The solution for y is given by y = -W_0(alpha) / alpha, where W_0 is the principal branch 
        # of the Lambert W function. For a positive real argument, W_0 is real and positive.
        w_val = lambertw(alpha, 0)
        y = -w_val.real / alpha
        
        z = y * 1j

        # Verification
        # LHS = z * i = iy * i = -y
        lhs = -y
        
        # RHS = i^z = i^(iy). We MUST use the logarithm branch for the corresponding k.
        log_i_k = 1j * (np.pi/2 + 2*k*np.pi)
        rhs = cmath.exp(z * log_i_k)

        print(f"For k = {k}:")
        print(f"Solution: z = {z:.8f}")
        print(f"The equation is: ({z:.6f}) * i = i^({z:.6f})")
        print(f"Verification: LHS = {lhs:.6f}, RHS = {rhs.real:.6f}. Match: {np.isclose(lhs, rhs.real)}\n")

solve_and_verify_complex_equation()