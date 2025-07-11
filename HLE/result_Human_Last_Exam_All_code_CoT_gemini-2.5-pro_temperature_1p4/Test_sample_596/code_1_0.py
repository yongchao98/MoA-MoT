# This script calculates the number of higher dimensional rooted forests of the
# Möbius band that do not simplicially collapse onto their root.

# This number, according to a theorem by Bernardi, de Mier, and Rué (2018),
# depends only on the homology of the underlying space, not its specific triangulation.
# The formula is N = |tilde_H_0(M, Z)|^2 + sum_{i>=1} (dim_Q H_i(M, Q))^2
# where M is the Möbius band.

# Step 1: Determine the size of the reduced 0-th homology group with integer coefficients.
# The Möbius band is path-connected, so its reduced 0-th homology group is trivial (0).
h0_tilde_abs = 0

# Step 2: Determine the dimensions of the higher homology groups with rational coefficients.
# The Möbius band is homotopy equivalent to a circle, so its first homology group H_1(M, Z) is Z.
# The dimension over the rationals Q is dim_Q(Z tensor Q) = dim_Q(Q) = 1.
dim_h1 = 1

# As a 2-manifold, the homology groups of the Möbius band are zero for dimensions 2 and higher.
# (Specifically, H_2 is zero because it is non-orientable).
# So, their dimensions are 0.
dim_h2 = 0
# Higher dimensions are also 0.

# Step 3: Apply the formula.
# The formula for a 2-dimensional surface simplifies to:
# N = |tilde_H_0(M, Z)|^2 + (dim_Q H_1(M, Q))^2 + (dim_Q H_2(M, Q))^2
result = h0_tilde_abs**2 + dim_h1**2 + dim_h2**2

# Print the final result along with the equation structure.
print("The number (N) of non-collapsing higher dimensional rooted forests is calculated using the homology of the Möbius band (M).")
print("The formula is: N = |tilde_H_0(M, Z)|^2 + (dim_Q H_1(M, Q))^2 + (dim_Q H_2(M, Q))^2")
print("\nPlugging in the values for the Möbius band:")
print(f"N = {h0_tilde_abs}^2 + {dim_h1}^2 + {dim_h2}^2")
print(f"N = {h0_tilde_abs**2} + {dim_h1**2} + {dim_h2**2}")
print(f"N = {result}")
