# This script calculates the properties of the three given CFGs and formats the output.

# --- Analysis of X1 ---
# X1 is the Hilbert scheme Hilb^11(A^3).
# Type: Scheme (S)
# Properties: separated (s)
# Dimension: For Hilb^d(A^n), the dimension is n * d.
n_1 = 3
d_1 = 11
dim_1 = n_1 * d_1
# Final profile for X1. It is reducible and not universally closed.
profile_1 = f"[S, s, {dim_1}]"

# --- Analysis of X2 ---
# X2 is the quotient stack [(A^4 \ V(xy-zw))/C*].
# Type: Deligne-Mumford stack (DM), as stabilizers are finite.
# Properties: separated (s), irreducible (irr).
# Dimension: For a quotient stack [U/G], the dimension is dim(U) - dim(G).
dim_space_2 = 4  # Dimension of A^4 \ V(xy-zw) is the dimension of A^4
dim_group_2 = 1  # Dimension of the group C*
dim_2 = dim_space_2 - dim_group_2
# Final profile for X2. It is not universally closed.
profile_2 = f"[DM, s, irr, {dim_2}]"

# --- Analysis of X3 ---
# X3 is the Picard stack Pic(C_0) for a genus 7 curve C_0.
# Type: Algebraic stack (A), as stabilizers are C*.
# Properties: separated (s).
# Dimension: The dimension of the Picard stack of a curve is the genus g.
g_3 = 7
dim_3 = g_3
# Final profile for X3. It is not irreducible and not universally closed.
profile_3 = f"[A, s, {dim_3}]"

# --- Combine and Print ---
# The final answer is the combination of the three profiles.
final_answer = f"{profile_1} {profile_2} {profile_3}"
print(final_answer)