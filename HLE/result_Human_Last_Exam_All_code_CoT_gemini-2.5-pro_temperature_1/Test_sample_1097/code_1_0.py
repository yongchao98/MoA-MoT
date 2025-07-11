import numpy as np
from scipy.optimize import fsolve

# The Banach Fixed-Point Theorem can be applied if we find a set M where our operator is a contraction.
# This leads to two conditions on the radius R of the set M = {u | ||u||_inf <= R}:
# 1. Contraction: exp(R) / 8 < 1  => R < ln(8)
# 2. Invariance (T(M) subset M): exp(R) / 8 <= R

# We can use Python to find the valid range for R.

# The upper bound for R comes from the contraction condition.
R_upper_bound = np.log(8)

# The lower bound comes from solving the equation exp(R) / 8 = R, or 8*R - exp(R) = 0.
# This equation has two roots. We are interested in the smaller positive root.
def g(R):
    """Function to find the root for."""
    return 8 * R - np.exp(R)

# We use fsolve to find the root. We need an initial guess.
# Let's guess a small positive number.
initial_guess = 0.1
R_lower_bound = fsolve(g, initial_guess)[0]

print("--- Analysis of the conditions on R ---")
print(f"1. Contraction condition requires R < ln(8), so R < {R_upper_bound:.4f}")
print(f"2. Invariance condition requires R >= {R_lower_bound:.4f}")
print(f"Therefore, a valid R must be in the interval [{R_lower_bound:.4f}, {R_upper_bound:.4f}).")

# We can choose any R in this interval. R=1 is a simple choice.
R_chosen = 1
is_valid = R_lower_bound <= R_chosen < R_upper_bound
print(f"\nLet's choose a simple value, R = {R_chosen}.")
print(f"Is R = {R_chosen} a valid choice? {is_valid}")

# Now we print the final answer: the definition of the set M.
# The prompt requires printing each number in the final equation.
boundary_val = 0
boundary_point_1 = 0
boundary_point_2 = 1

print("\n--- Final Answer ---")
print("The set M you should define is:")
# The following print statement constructs the definition of M as requested.
print(f"M = {{u in C[0, 1] | u({boundary_point_1}) = u({boundary_point_2}) = {boundary_val} and max|u(x)| <= {R_chosen}}}")
print("where C[0, 1] is the space of continuous functions on [0, 1] and max|u(x)| is the supremum norm.")
