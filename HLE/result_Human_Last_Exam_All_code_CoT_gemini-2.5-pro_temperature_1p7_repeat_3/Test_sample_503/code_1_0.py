import math

# The problem is to find the number of internal adjunctions from [m] to [n].
# Here m = 23 and n = 37.
m = 23
n = 37

# As explained, this number is equivalent to the number of non-decreasing maps
# L: [m] -> [n] such that L(0) = 0.
# This corresponds to choosing k=m values (for L(1),...,L(m)) from a set
# of s=n+1 possible values ({0, ..., n}) with replacement.
# The formula is C(s + k - 1, k) = C((n+1) + m - 1, m) = C(n+m, m).
comb_n = n + m
comb_k = m

# Calculate the binomial coefficient C(60, 23).
result = math.comb(comb_n, comb_k)

# Print the final result in a descriptive way.
print(f"The number of internal adjunctions from [{m}] to [{n}] is given by the binomial coefficient C({comb_n}, {comb_k}).")
print(f"C({comb_n}, {comb_k}) = {result}")
