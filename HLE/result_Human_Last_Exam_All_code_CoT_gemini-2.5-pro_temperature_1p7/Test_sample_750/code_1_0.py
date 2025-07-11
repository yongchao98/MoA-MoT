import numpy as np
from scipy.optimize import root_scalar

def solve_for_y(k):
    """
    Numerically solves the transcendental equation -y = exp(-y*(pi/2 + 2*k*pi)) for y.
    This equation has real solutions for y only for k < 0.
    """
    if k >= 0:
        return None
    
    # The constant in the exponent
    beta_k = np.pi/2 + 2*k*np.pi
    
    # The function whose root we want to find: f(y) = y + exp(-y * beta_k)
    func_to_solve = lambda y: y + np.exp(-y * beta_k)
    
    # For k < 0, beta_k is negative. We can prove that a root exists between -1 and 0.
    # At y=0, func(0) = 1.
    # At y=-1, func(-1) = -1 + exp(beta_k) < 0 since beta_k < 0.
    # So we can safely bracket the root in [-1, 0].
    solution = root_scalar(func_to_solve, bracket=[-1, 0])
    
    return solution.root

# Main script execution
print("Solving the complex equation: z * i = i^z")
print("The solutions for z are:\n")

# --- Real Solutions ---
print("1. Two real solutions exist:")
# Solution z=1
z1 = 1
print(f"   - For z = {z1}:")
print(f"     The equation is {z1} * i = i^{z1}\n")

# Solution z=-1
z2 = -1
print(f"   - For z = {z2}:")
print(f"     The equation is {z2} * i = i^({z2})\n")

# --- Purely Imaginary Solutions ---
print("2. An infinite series of purely imaginary solutions (z = iy) also exists.")
print("   These correspond to different branches of the complex logarithm (for integers k < 0).\n")

# Calculate and print the first 3 imaginary solutions
for k in range(-1, -4, -1):
    y_k = solve_for_y(k)
    
    # Format the imaginary number for printing
    z_k_str = f"({y_k:.6f}i)"
    
    print(f"   - For k = {k}, a solution is z ≈ {y_k:.6f}i:")
    print(f"     The equation is {z_k_str} * i = i^{z_k_str}\n")
