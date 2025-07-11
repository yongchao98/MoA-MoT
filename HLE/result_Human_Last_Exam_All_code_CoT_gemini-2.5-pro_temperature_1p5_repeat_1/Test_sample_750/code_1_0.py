import numpy as np
from scipy.special import lambertw

def solve_and_verify():
    """
    Finds and verifies solutions for z*i = i^z.
    """
    print("Solving for z in the complex equation z*i = i^z\n")
    
    # Case 1: Real solutions
    print("Found two real solutions: z = 1 and z = -1.\n")
    print("Verification for z = 1:")
    z1 = 1
    lhs1 = z1 * 1j
    # For z=1, i^1 = i, which corresponds to any branch k of log(i).
    # For instance, principal branch k=0: exp(1j*pi/2) = i
    rhs1 = 1j 
    print(f"z = {z1}")
    print(f"LHS = z * i = {z1} * {1j} = {lhs1}")
    print(f"RHS = i^z = i^1 = {rhs1}")
    print(f"Equation: {lhs1} = {rhs1}\n")

    print("Verification for z = -1:")
    z2 = -1
    lhs2 = z2 * 1j
    # For z=-1, i^-1 = 1/i = -i, corresponds to any branch k.
    # For instance, principal branch k=0: exp(-1j*pi/2) = -i
    rhs2 = -1j
    print(f"z = {z2}")
    print(f"LHS = z * i = {z2} * {1j} = {lhs2}")
    print(f"RHS = i^z = i^-1 = {rhs2}")
    print(f"Equation: {lhs2} = {rhs2}\n")

    # Case 2: Purely imaginary solutions
    print("Found a series of purely imaginary solutions z = yi for k <= -1.")
    print("These solutions are of the form z = -a*i, where a = W(D)/D and D = (4m-1)pi/2 for m >= 1.\n")
    
    for m in range(1, 4):
        k = -m
        # D_m = -(pi/2 + 2*k*pi) = (4m-1)pi/2
        D = (4 * m - 1) * np.pi / 2
        
        # Solve for a using the Lambert W function
        # The result of lambertw is complex, so we take the real part.
        a = np.real(lambertw(D) / D)
        
        z = -1j * a
        
        # Verification
        lhs = z * 1j
        
        # For the RHS, we must use the corresponding branch of log(i)
        # C_k = pi/2 + 2*k*pi
        C_k = np.pi/2 + 2*k*np.pi
        rhs = np.exp(z * 1j * C_k)

        print(f"For m = {m} (k = {k}):")
        print(f"z = {z}")
        print(f"LHS = z * i = ({z}) * {1j} = {lhs}")
        print(f"RHS = i^z (using k={k} branch) = exp(z * i * ({C_k:.4f})) = {rhs}")
        # The final equation demonstrates the equality of the numerical results
        print(f"Final Equation: {np.real(lhs):.12f} = {np.real(rhs):.12f}\n")

solve_and_verify()