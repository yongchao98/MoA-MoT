# This script determines the properties of three categories fibered in groupoids
# and prints the result in the specified format.

# Analysis of X1
# X1 is the Hilbert scheme of 11 points in A^3.
# Type: S (Scheme)
# Separated: s
# Not universally closed, not irreducible.
# Dimension is based on the principal component: dim(A^3) * degree
dim_X1 = 3 * 11
profile_X1 = f"[S, s, {dim_X1}]"

# Analysis of X2
# X2 is the quotient stack [(A^4 \ V(xy-zw))/C*].
# Type: DM (Deligne-Mumford stack due to finite stabilizers)
# Separated: s
# Irreducible: irr
# Not universally closed.
# Dimension is dim(space) - dim(group)
dim_X2 = 4 - 1
profile_X2 = f"[DM, s, irr, {dim_X2}]"

# Analysis of X3
# X3 is the Picard stack of a genus 7 curve.
# Type: A (Algebraic stack due to Gm stabilizers)
# Separated: s
# Not universally closed, not irreducible.
# Dimension is genus + dim(stabilizer)
dim_X3 = 7 + 1
profile_X3 = f"[A, s, {dim_X3}]"

# Combine the profiles into the final answer string
final_answer = f"{profile_X1} {profile_X2} {profile_X3}"

print(final_answer)