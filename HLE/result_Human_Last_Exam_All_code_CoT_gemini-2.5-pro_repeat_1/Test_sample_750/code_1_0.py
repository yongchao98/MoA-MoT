import numpy as np
import cmath
from scipy.special import lambertw

def solve_and_verify():
    """
    Solves the complex equation z*i = i^z and verifies the solutions.
    The final equation for each solution is printed by showing that LHS = RHS.
    """
    print("The solutions to the equation z * i = i^z are found to be:")
    
    # --- Real Solutions ---
    print("\n1. Real solutions:")
    real_solutions = [1.0, -1.0]
    for z_val in real_solutions:
        z = complex(z_val)
        lhs = z * 1j
        
        # For verification, we need to choose the correct branch of log(i).
        # For z=1, i = i^1 = exp(1*log(i)). arg(i) = arg(exp(log(i))).
        # pi/2 = (pi/2 + 2k*pi) mod 2pi. k=0 works.
        # For z=-1, -i = i^-1 = exp(-1*log(i)). arg(-i) = arg(exp(-log(i))).
        # -pi/2 = -(pi/2 + 2k*pi) mod 2pi. -pi/2 = -pi/2 - 2k*pi. Any k works.
        # We use the principal branch (k=0) for simplicity in verification where possible.
        log_i_principal = 1j * np.pi / 2
        rhs = cmath.exp(z * log_i_principal)
        # For z=-1, this gives i, not -i. We need another branch or use i**z directly.
        # Python's ** operator for complex numbers uses the principal branch formula: a**b = exp(b*log(a))
        rhs = (1j)**z

        print(f"\nFor z = {z.real}:")
        print(f"  LHS: z * i = ({z.real}) * i = {lhs.real:.4f}{lhs.imag:+.4f}j")
        print(f"  RHS: i^z = i^({z.real}) = {rhs.real:.4f}{rhs.imag:+.4f}j")
        if cmath.isclose(lhs, rhs):
            print("  Verification: LHS = RHS")
        else:
            print("  Verification: LHS != RHS (check log branch)")


    # --- Purely Imaginary Solutions ---
    print("\n2. An infinite series of purely imaginary solutions:")
    print("   z_j = -i * W_0(c_j) / c_j, where c_j = (2j - 0.5)*pi for j=1, 2, ...")
    print("   W_0 is the principal branch of the Lambert W function.")
    
    num_solutions_to_show = 5
    print(f"\n   The first {num_solutions_to_show} imaginary solutions are:")
    
    for j in range(1, num_solutions_to_show + 1):
        c_j = (2 * j - 0.5) * np.pi
        w_val = lambertw(c_j, k=0) # W_0(c_j)
        z = -1j * w_val / c_j

        # Verification
        lhs = z * 1j
        # For RHS, we must use the correct branch of log(i) that gave this solution.
        # These solutions correspond to k = -j in the log definition log(i) = i(pi/2 + 2k*pi)
        k = -j
        log_i_branch = 1j * (np.pi/2 + 2 * k * np.pi)
        rhs = cmath.exp(z * log_i_branch)
        
        print(f"\nFor j = {j} (c_{j} = {c_j:.4f}):")
        print(f"  z_{j} = {z.real:.4f}{z.imag:+.4f}j")
        print(f"  LHS: z_{j} * i = {lhs.real:.4f}{lhs.imag:+.4f}j")
        print(f"  RHS: i^(z_{j}) = {rhs.real:.4f}{rhs.imag:+.4f}j  (using log branch k={k})")
        if cmath.isclose(lhs, rhs):
            print("  Verification: LHS = RHS")
        else:
            print("  Verification: LHS != RHS")

if __name__ == '__main__':
    solve_and_verify()