import collections

# Based on the reasoning above, the sequence D_g can be determined.
# D_g is the minimal degree of a cover of the moduli stack A_g
# that trivializes the obstruction to finding a symmetric line bundle
# representing the universal polarization.

# The obstruction lives in H^1(Sp(2g, Z), (Z/2Z)^{2g}).
# For g=1, this group is 0.
d1 = 1
# For g=2, this group is 0.
d2 = 1
# For g=3, this group is Z/2Z.
d3 = 2
# For g=4, this group is Z/2Z.
d4 = 2

# The sequence is formed by these values.
sequence = [d1, d2, d3, d4]

print("The sequence of the first 4 terms of D_g is derived from group cohomology.")
print(f"For g = 1, the obstruction group is trivial, so no cover is needed. D_1 = {sequence[0]}")
print(f"For g = 2, the obstruction group is trivial, so no cover is needed. D_2 = {sequence[1]}")
print(f"For g = 3, the obstruction group is Z/2Z, requiring a degree 2 cover. D_3 = {sequence[2]}")
print(f"For g = 4, the obstruction group is Z/2Z, requiring a degree 2 cover. D_4 = {sequence[3]}")

# The problem asks to write the sequence of the first 4 terms.
# The following print statement provides the final answer in a clear format.
print("\nThe sequence D_g for g=1, 2, 3, 4 is:")
print(", ".join(map(str, sequence)))