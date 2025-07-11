import numpy as np
from scipy.special import lambertw

def solve_and_verify():
    """
    Solves the equation z*i = i^z for real and purely imaginary solutions
    and verifies them.
    """
    # Define i as a complex number
    i = 1j

    print("Solving for z in z*i = i^z\n")

    # --- Real Solutions (z = x) ---
    print("--- Real Solutions ---")
    
    # Solution z = 1
    z1 = 1
    lhs1 = z1 * i
    rhs1 = i**z1
    print(f"Solution: z = {z1}")
    # Final equation format as requested
    print(f"Verification: The equation is {z1} * i = i^{z1}")
    print(f"LHS = {z1} * {i} = {lhs1}")
    print(f"RHS = {i}^{z1} = {rhs1}")
    print(f"Result is {'True' if np.isclose(lhs1, rhs1) else 'False'}\n")

    # Solution z = -1
    z2 = -1
    lhs2 = z2 * i
    rhs2 = i**z2
    print(f"Solution: z = {z2}")
    # Final equation format as requested
    print(f"Verification: The equation is {z2} * i = i^{z2}")
    print(f"LHS = {z2} * {i} = {lhs2}")
    print(f"RHS = {i}^{z2} = {rhs2}")
    print(f"Result is {'True' if np.isclose(lhs1, rhs1) else 'False'}\n")

    # --- Purely Imaginary Solutions (z = yi) ---
    # The equation for y is -y = exp(-y * alpha_k) where alpha_k = pi/2 + 2k*pi
    # This only has solutions for k <= -1.
    # The solution for y is given by y = -W(-alpha_k) / (-alpha_k) where W is the Lambert W function.
    # So z = yi = -i * W(-alpha_k) / (-alpha_k)
    
    print("--- Purely Imaginary Solutions ---")
    print("This family of solutions is z_k = yi, where y solves -y = exp(-y*(pi/2 + 2k*pi)) for k <= -1.")
    print("Calculating the first 3 solutions (for k = -1, -2, -3):\n")

    for k in range(-1, -4, -1):
        alpha_k = np.pi/2 + 2*k*np.pi
        
        # We solve t = exp(t*alpha_k) for t = -y.
        # Let alpha_k = -beta. t = exp(-t*beta). t*beta*exp(t*beta) = beta.
        # t*beta = W(beta), so t = W(beta)/beta.
        # beta = -alpha_k
        beta = -alpha_k
        
        # The argument to lambertw is positive, so the result is real.
        t = np.real(lambertw(beta)) / beta
        y = -t
        z = y * i
        
        print(f"Solution for k = {k}:")
        print(f"z = {z:.8f}")

        # Verification
        lhs = z * i
        
        # For the RHS, we must use the correct branch of the logarithm
        # i^z = exp(z * log(i)) = exp(z * i * alpha_k)
        rhs = np.exp(z * i * alpha_k)

        # Final equation format as requested
        print(f"Verification: The equation is ({z:.4f}) * i = i^({z:.4f})")
        print(f"LHS = z * i = {lhs:.8f}")
        print(f"RHS = exp(z * i * (pi/2 + 2*({k})*pi)) = {rhs:.8f}")
        print(f"Result is {'True' if np.isclose(lhs, rhs) else 'False'}\n")

solve_and_verify()