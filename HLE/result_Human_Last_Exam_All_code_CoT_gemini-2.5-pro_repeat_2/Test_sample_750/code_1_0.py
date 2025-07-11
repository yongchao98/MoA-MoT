import numpy as np
from scipy.special import lambertw

def solve_z(k, j):
    """
    Calculates a solution z for the equation z*i = i^z based on integer indices k and j.
    k: integer for the branch of the logarithm of i.
    j: integer for the branch of the Lambert W function.
    """
    # Calculate C_k = (pi/2 + 2*k*pi)
    C_k = np.pi / 2 + 2 * k * np.pi
    
    # Argument for the Lambert W function
    arg_W = -C_k
    
    # Calculate a = W_j(-C_k) / C_k
    try:
        a = lambertw(arg_W, k=j) / C_k
    except (ValueError, TypeError):
        return None # In case of invalid inputs for lambertw
    
    # Calculate z = a * i
    z = a * 1j
    
    return z

def main():
    """
    Calculates and prints various solutions to z*i = i^z.
    """
    print("The equation is z * i = i^z.")
    print("Solutions are of the form z = i * W_j(-(pi/2 + 2*k*pi)) / (pi/2 + 2*k*pi)")
    print("where k is the integer index for the logarithm branch and j is for the Lambert W function branch.\n")

    print("Two special real solutions are:")
    # We can verify z=1 and z=-1 are solutions
    print(f"z = {1.0:.1f}")
    print(f"z = {-1.0:.1f}\n")

    print("Examples of purely imaginary solutions (for k <= -1, j = 0):")
    for k in range(-1, -4, -1):
        z = solve_z(k, 0)
        print(f"For k={k}, j=0:  z = {z:.6f}")
    
    print("\nExamples of complex solutions (for k >= 0):")
    # For k=0, let's check a few branches j
    k = 0
    for j in range(-2, 3):
        z = solve_z(k, j)
        print(f"For k={k}, j={j}: z = {z.real:.6f} + {z.imag:.6f}i")

    # For k=1, let's check a few branches j
    k = 1
    print()
    for j in range(-2, 3):
        z = solve_z(k, j)
        print(f"For k={k}, j={j}: z = {z.real:.6f} + {z.imag:.6f}i")

if __name__ == "__main__":
    main()