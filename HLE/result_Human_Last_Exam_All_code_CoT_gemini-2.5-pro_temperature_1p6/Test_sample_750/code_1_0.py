import numpy as np
from scipy.optimize import fsolve
import warnings

# Suppress warning from fsolve if it doesn't converge, though it's expected to here.
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

def solve_and_print():
    """
    Finds and prints solutions to the equation z*i = i^z.
    """
    print("Solving the complex equation: z * i = i^z")
    print("="*40)

    # Case 1: Real solutions (z = x + 0i)
    # The analysis shows two real solutions.
    print("Found two real solutions:")
    z1 = 1.0
    z2 = -1.0
    print(f"z = {z1}")
    print(f"z = {z2}")
    
    # Verification for real solutions
    # For z = 1: 1*i = i; i^1 = i. Correct.
    # For z = -1: -1*i = -i; i^-1 = 1/i = -i. Correct.
    
    print("-" * 40)
    
    # Case 2: Purely imaginary solutions (z = 0 + yi)
    # This leads to the equation: -y = exp(-y * A_k)
    # where A_k = pi/2 + 2*k*pi.
    # Solutions exist for k <= -1.
    print("Found an infinite series of purely imaginary solutions z = iy.")
    print("The values y are roots of -y = exp(-y * (pi/2 + 2*k*pi)) for k <= -1.")
    print("Here are the solutions for the first 5 values of k:")
    
    # Loop for k from -1 to -5
    for k in range(-1, -6, -1):
        # A_k is the constant for a given k
        A_k = np.pi / 2 + 2 * k * np.pi
        
        # We need to find the root of the function f(y) = y + exp(-y * A_k)
        def func_to_solve(y):
            return y + np.exp(-y * A_k)

        # Analysis shows the root y_k is in (-1, 0). 
        # A good initial guess helps the solver.
        initial_guess = -0.1
        y_solution = fsolve(func_to_solve, initial_guess)[0]
        
        # The solution is z = iy
        z_imag = y_solution * 1j
        
        print(f"For k = {k}, a solution is z = {z_imag:.8f}")

    print("="*40)
    print("Analysis shows that z=1 and z=-1 are the only solutions on the unit circle.")
    print("Other solutions with non-zero real and imaginary parts may exist but are not found by this analysis.")

solve_and_print()