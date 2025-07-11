# For X_1, the Hilbert scheme of 11 points in A^3.
# Its dimension is d * n, where d is the number of points and n is the dimension of the ambient space.
d_1 = 11
n_1 = 3
dim_1 = d_1 * n_1
# Properties: Scheme, separated, not irreducible, not universally closed.
profile_1 = f"[S,s,{dim_1}]"

# For X_2, the quotient stack [(\mathbb{A}^4 \setminus V(xy-zw))/\mathbb{C}^*].
# Its dimension is dim(total space) - dim(group).
dim_space_2 = 4
dim_group_2 = 1
dim_2 = dim_space_2 - dim_group_2
# Properties: Scheme, separated, irreducible, not universally closed.
profile_2 = f"[S,s,irr,{dim_2}]"

# For X_3, the Picard stack of a genus 7 curve.
# Its dimension is g - 1, where g is the genus.
g_3 = 7
dim_stab_3 = 1
dim_3 = g_3 - dim_stab_3
# Properties: Algebraic stack, separated, not irreducible, not universally closed.
profile_3 = f"[A,s,{dim_3}]"

# Combine the profiles into the final answer string.
final_answer = f"{profile_1} {profile_2} {profile_3}"

print(final_answer)