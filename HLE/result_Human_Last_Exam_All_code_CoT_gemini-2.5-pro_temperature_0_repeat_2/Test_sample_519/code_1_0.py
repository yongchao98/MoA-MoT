# For X1, the Hilbert scheme of 11 points in A^3.
# The dimension is n * d.
dim1 = 3 * 11
profile1 = f"[S,s,{dim1}]"

# For X2, the quotient space [(A^4 \ V(xy-zw))/C*].
# The dimension is dim(total space) - dim(group).
dim2 = 4 - 1
profile2 = f"[S,s,irr,{dim2}]"

# For X3, the Picard scheme of a genus 7 curve.
# The dimension is the genus of the curve.
dim3 = 7
profile3 = f"[S,s,{dim3}]"

# Combine the profiles into a single string.
final_answer = f"{profile1} {profile2} {profile3}"

print(final_answer)