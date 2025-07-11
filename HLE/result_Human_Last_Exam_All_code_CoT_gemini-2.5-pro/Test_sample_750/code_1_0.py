import numpy as np
from scipy.special import lambertw

def solve_and_verify():
    """
    Finds and verifies solutions to the complex equation z*i = i^z.
    """
    print("Solving for z in the equation z*i = i^z\n")

    # --- Real Solutions ---
    print("--- Real Solutions ---")
    
    # Solution z = 1
    z1 = 1 + 0j
    # For z=1, i^1 is typically evaluated using the principal branch k=0
    # alpha = pi/2 + 2*k*pi. For k=0, alpha=pi/2.
    lhs1 = z1 * 1j
    rhs1 = 1j**z1 # Python's ** uses the principal branch k=0
    print(f"Solution 1: z = {z1}")
    print(f"Verification for z=1:")
    print(f"  LHS = z * i = {z1} * i = {lhs1}")
    print(f"  RHS = i^z = i^1 = {rhs1}")
    print(f"  Result: {'MATCH' if np.isclose(lhs1, rhs1) else 'NO MATCH'}\n")

    # Solution z = -1
    z2 = -1 + 0j
    # For z=-1, i^-1 is typically evaluated using the principal branch k=0
    lhs2 = z2 * 1j
    rhs2 = 1j**z2 # Python's ** uses the principal branch k=0
    print(f"Solution 2: z = {z2}")
    print(f"Verification for z=-1:")
    print(f"  LHS = z * i = {z2} * i = {lhs2}")
    print(f"  RHS = i^z = i^-1 = {rhs2}")
    print(f"  Result: {'MATCH' if np.isclose(lhs2, rhs2) else 'NO MATCH'}\n")
    
    # --- Purely Imaginary Solutions ---
    print("--- Purely Imaginary Solutions (z = iy) ---")
    print("These solutions z_k = i*y_k exist for integer k <= -1.")
    print("y_k is the solution to -y = exp(-y * alpha), where alpha = pi/2 + 2*k*pi.")
    print("We calculate and verify the first 5 such solutions.\n")

    for k in range(-1, -6, -1):
        # alpha corresponds to the branch of the logarithm used for i^z
        alpha = np.pi/2 + 2*k*np.pi
        
        # We solve -y = exp(-y*alpha) using the Lambert W function.
        # The equation is equivalent to y * exp(y*alpha) = -1.
        # Let u = y*alpha. Then (u/alpha)*exp(u) = -1 => u*exp(u) = -alpha.
        # So, u = W(-alpha), which means y*alpha = W(-alpha).
        # y = W(-alpha) / alpha
        
        # For k <= -1, alpha is negative, so -alpha is positive.
        # lambertw(-alpha) will return a real number.
        y_k = (lambertw(-alpha)).real / alpha
        z_k = 0 + y_k*1j
        
        # Verification
        # LHS = z_k * i = (i*y_k)*i = -y_k
        lhs_k = -y_k
        # RHS = (i^z)_k = exp(z_k * log_k(i)) = exp((i*y_k) * (i*alpha)) = exp(-y_k*alpha)
        rhs_k = np.exp(-y_k * alpha)
        
        print(f"Solution for k={k}:")
        print(f"  z_{k} = {z_k:.6f}")
        print(f"  Verification for z_{k}:")
        print(f"    LHS (-y_k) = {lhs_k:.6f}")
        print(f"    RHS (exp(-y_k*alpha)) = {rhs_k:.6f}")
        print(f"    Result: {'MATCH' if np.isclose(lhs_k, rhs_k) else 'NO MATCH'}\n")

solve_and_verify()