import numpy as np
from scipy.special import lambertw

def solve_and_print_zi_iz():
    """
    Finds and prints the solutions to the equation z*i = i^z.
    The solutions are found by representing z as x+iy and solving the resulting
    system of equations for real and imaginary parts. This leads to two real solutions
    and an infinite series of purely imaginary solutions which can be expressed
    using the Lambert W function.
    """
    print("The equation to solve is z*i = i^z for all complex numbers z.")
    print("\n--- Real Solutions ---")
    print("There are two real solutions:")
    print("z = 1")
    print("z = -1")

    print("\n--- Purely Imaginary Solutions ---")
    print("There is an infinite series of purely imaginary solutions of the form:")
    print("z_m = -i * W(2*pi*m - pi/2) / (2*pi*m - pi/2) for m = 1, 2, 3, ...")
    print("where W is the principal branch of the Lambert W function.")
    print("\nHere are the first 5 imaginary solutions:")

    for m in range(1, 6):
        # Argument for the Lambert W function
        d_m = 2 * np.pi * m - np.pi / 2
        
        # The lambertw function returns a complex number; for a positive real argument,
        # the principal branch W_0 is real.
        w_val = np.real(lambertw(d_m))
        
        # Calculate y-part of the solution
        y_m = -w_val / d_m
        
        # The solution z_m is purely imaginary
        z_m = y_m * 1j

        print(f"\nFor m={m}:")
        print(f"  The argument of W is d_{m} = 2*pi*{m} - pi/2 = {d_m:.4f}")
        print(f"  The value of W(d_{m}) is {w_val:.4f}")
        # The final equation output as requested
        print(f"  The final equation for z_{m} is -i * W({d_m:.4f}) / {d_m:.4f}")
        print(f"  The solution is z_{m} = {z_m.real:.4f} {z_m.imag:+.4f}i")

solve_and_print_zi_iz()