# Properties for X1: Hilb^11(A^3)
dim1 = 11 * 3
# It is a Scheme (S), separated (s), but not universally closed or irreducible.
# Dimension is 33.
x1_profile = f"[S,s,{dim1}]"

# Properties for X2: [(A^4 \ V(xy-zw))/C*]
dim2 = 4 - 1
# It is a DM-stack (DM), separated (s), and irreducible (irr).
# It is not universally closed.
# Dimension is 3.
x2_profile = f"[DM,s,irr,{dim2}]"

# Properties for X3: Pic_C0 for g=7 curve C0
g = 7
dim_gm = 1
dim3 = g + dim_gm
# It is an Algebraic stack (A), separated (s), but not universally closed or irreducible.
# Dimension is 8.
x3_profile = f"[A,s,{dim3}]"

# Combine the profiles into a single string
final_output = f"{x1_profile} {x2_profile} {x3_profile}"

print(final_output)