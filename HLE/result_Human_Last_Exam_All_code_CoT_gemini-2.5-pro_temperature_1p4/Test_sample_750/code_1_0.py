import numpy as np
from scipy.special import lambertw

def solve_z_i_eq_iz():
    """
    Solves for z in the complex equation z*i = i^z.
    """
    print("The equation z*i = i^z has two real solutions and an infinite series of purely imaginary solutions.")
    
    # Real solutions
    print("\nReal solutions:")
    z1 = 1.0
    z2 = -1.0
    print(f"z = {z1} + 0.0i")
    print(f"z = {z2} + 0.0i")
    
    # Purely imaginary solutions
    print("\nPurely imaginary solutions of the form z = iy:")
    print("The solutions are of the form z_m = i * y_m for m = 1, 2, 3, ...")
    print("where y_m = -W((4m-1)pi/2) / ((4m-1)pi/2)")
    print("W is the principal branch of the Lambert W function.")
    print("\nFirst 5 purely imaginary solutions:")
    
    for m in range(1, 6):
        # B_m = (4m-1)pi/2
        # We pass it as a complex number to ensure numerical stability in the library functions.
        B_m = complex((4 * m - 1) * np.pi / 2, 0)
        
        # Calculate the Lambert W function for B_m
        W_Bm = lambertw(B_m, k=0)
        
        # Calculate y_m
        y_m = -W_Bm / B_m
        
        # The solution is z_m = i * y_m
        z_m = 1j * y_m
        
        # Outputting the numbers in the final form
        # The imaginary part of z_m will be the real part of y_m
        print(f"For m={m}, z_{m} = {z_m.real:.6f} + {z_m.imag:.6f}i")

if __name__ == '__main__':
    solve_z_i_eq_iz()