import numpy as np
from scipy.special import lambertw

def solve_and_verify():
    """
    Finds and verifies solutions to the equation z*i = i^z.
    """
    print("Solving the complex equation: z * i = i^z")
    print("="*40)

    # --- Real Solutions ---
    print("Found two real solutions:")
    
    # Solution z = 1
    z1 = 1.0
    lhs1 = z1 * 1j
    # For z=1 or z=-1, the principal value of i^z is sufficient
    rhs1 = 1j**z1
    print(f"\n1. Solution z = {z1}")
    print(f"   Left side (z*i): {lhs1}")
    print(f"   Right side (i^z): {rhs1}")
    print(f"   Verification (LHS == RHS): {np.isclose(lhs1, rhs1)}")

    # Solution z = -1
    z2 = -1.0
    lhs2 = z2 * 1j
    rhs2 = 1j**z2
    print(f"\n2. Solution z = {z2}")
    print(f"   Left side (z*i): {lhs2}")
    print(f"   Right side (i^z): {rhs2}")
    print(f"   Verification (LHS == RHS): {np.isclose(lhs2, rhs2)}")
    
    print("\n" + "="*40)

    # --- Purely Imaginary Solutions ---
    print("Found an infinite series of purely imaginary solutions z_k = y_k * i")
    print("for negative integers k = -1, -2, -3, ...")
    print("The formula for y_k is: y_k = W(-theta_k) / theta_k")
    print("where theta_k = pi/2 + 2*k*pi, and W is the Lambert W function.")
    print("\nCalculating the first 5 imaginary solutions:")

    for k in range(-1, -6, -1):
        theta_k = np.pi/2 + 2*k*np.pi
        
        # Calculate y_k using the principal branch of the Lambert W function.
        # The argument -theta_k is positive, so the result is a single real number.
        y_k = lambertw(-theta_k).real / theta_k
        z_k = y_k * 1j
        
        # Verification
        lhs_k = z_k * 1j
        
        # For the RHS, we must use the specific branch of log(i) corresponding to k.
        # i^z = exp(z * log(i)), and log(i) = i * theta_k
        log_i_k = 1j * theta_k
        rhs_k = np.exp(z_k * log_i_k)

        print(f"\nFor k={k}:")
        print(f"   Solution z_{k} = {z_k:.8f}")
        print(f"   Left side (z*i): {lhs_k:.8f}")
        print(f"   Right side (i^z using log branch k={k}): {rhs_k:.8f}")
        print(f"   Verification (LHS == RHS): {np.isclose(lhs_k, rhs_k)}")

if __name__ == '__main__':
    solve_and_verify()
