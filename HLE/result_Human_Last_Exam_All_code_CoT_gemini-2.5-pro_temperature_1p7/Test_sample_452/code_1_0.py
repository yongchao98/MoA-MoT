import math

# The problem asks for the constant 'b' in the asymptotic formula C(n) ~ b * n^(3/2),
# where C(n) is the expected cover-and-return time for a random walk on a random tree.

# Step 1: The expected cover-and-return time from a vertex v in a tree is given by
# the identity: E_v[T_CR] = sum_{u in V} d(u) * dist(v, u),
# where d(u) is the degree of vertex u and dist(v, u) is the distance between v and u.

# Step 2: C(n) is the expectation of this quantity, averaged over all random trees T
# on n vertices and all starting vertices v.
# C(n) = E_T [ (1/n) * sum_{v in V} sum_{u in V} d(u) * dist(v, u) ]
# The double summation is the degree-distance DD(T) of the tree.
# C(n) = (1/n) * E_T[DD(T)]

# Step 3: A known result for the expected degree-distance of a uniform random tree is:
# E_T[DD(T)] ~ sqrt(2 * pi) * n^(5/2)

# Step 4: Substituting this into the expression for C(n) gives:
# C(n) ~ (1/n) * (sqrt(2 * pi) * n^(5/2)) = sqrt(2 * pi) * n^(3/2)

# Step 5: By comparing this with C(n) ~ b * n^(3/2), we find that b = sqrt(2 * pi).

# The Python code below calculates this value.
number_2 = 2
pi = math.pi
b = math.sqrt(number_2 * pi)

print(f"The exact value of the constant b is the square root of (2 * pi).")
print(f"Equation: b = sqrt({number_2} * {pi})")
print(f"Numerical value: b = {b}")
