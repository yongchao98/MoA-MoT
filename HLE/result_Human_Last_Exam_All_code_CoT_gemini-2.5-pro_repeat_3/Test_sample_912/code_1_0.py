import sympy as sp

# Define the symbols
mu, mu_0, N, w, g, x1, x2, I1, I2 = sp.symbols('mu mu_0 N w g x_1 x_2 I_1 I_2')

# Construct the expression for the work done
# This is based on the derivation W = -1/2 * (L(x2) - L(x1)) * (I2**2 - I1**2)
# where dL/dx = N**2 * w * (mu - mu_0) / g
# so L(x2) - L(x1) = (dL/dx) * (x2 - x1)
dL_dx = (N**2 * w * (mu - mu_0)) / g
L_diff = dL_dx * (x2 - x1)
I_sq_diff = I2**2 - I1**2
Work = -sp.Rational(1, 2) * L_diff * I_sq_diff

# Rearrange to match the format in option D
Work_final = -((mu - mu_0) / (2 * g)) * N**2 * w * (x2 - x1) * (I2**2 - I1**2)

# Print the final equation
# The prompt requests to "output each number in the final equation"
# The following print statement shows the structure and components of the final formula.
print("The work done by the current source for each cycle is W:")
print(f"W = - ( (mu - mu_0) / (2*g) ) * N^2 * w * (x_2 - x_1) * (I_2^2 - I_1^2)")

# To be more explicit about each component as requested:
print("\nBreaking down the equation:")
print(f"Sign: -")
print(f"Permeability Term: (mu - mu_0)")
print(f"Gap Geometry Term: 1 / (2*g)")
print(f"Coil Term: N^2")
print(f"Area Term: w")
print(f"Displacement Term: (x_2 - x_1)")
print(f"Current Term: (I_2^2 - I_1^2)")
