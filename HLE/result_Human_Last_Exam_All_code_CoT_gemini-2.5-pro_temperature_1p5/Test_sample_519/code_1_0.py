# This script determines the properties of the three given categories fibered in groupoids
# and prints their profiles in the specified format.

# --- Case X1: Hilbert scheme ---
# Type: Scheme (S)
# Properties: separated (s)
# Dimension: For Hilb^d(A^n), the dimension is n * d.
n1 = 3
d1 = 11
dim1 = n1 * 11
profile1 = f"[S, s, {dim1}]"

# --- Case X2: Quotient Stack ---
# Type: Deligne-Mumford stack (DM)
# Properties: separated (s), irreducible (irr)
# Dimension: For [U/G], the dimension is dim(U) - dim(G).
dim_U = 4
dim_G = 1
dim2 = dim_U - dim_G
profile2 = f"[DM, s, irr, {dim2}]"

# --- Case X3: Picard Stack ---
# Type: Algebraic stack (A)
# Properties: separated (s)
# Dimension: The dimension of the Picard stack of a curve is the genus of the curve.
genus = 7
dim3 = genus
profile3 = f"[A, s, {dim3}]"

# Combine and print the final result
final_answer = f"{profile1} {profile2} {profile3}"
print(final_answer)