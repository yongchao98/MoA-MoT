# Based on the derivation, the problem of maximizing H(x,y,z,s1,s2)
# reduces to maximizing H(y,z).

# The maximum value of H(y,z) is bounded by the sum of individual entropies:
# H(y,z) <= H(y) + H(z)

# The constraints give us the maximum values for H(y) and H(z).
H_y = 1
H_z = 1

# The upper bound for the total entropy is the sum of these maximums.
# A construction has been shown to exist that achieves this bound.
max_H = H_y + H_z

# The final equation is derived from the simplification and bound.
# H_max = H(y) + H(z)
print("The maximal entropy can be calculated as follows:")
print(f"H_max = H(y) + H(z) = {H_y} + {H_z} = {max_H}")