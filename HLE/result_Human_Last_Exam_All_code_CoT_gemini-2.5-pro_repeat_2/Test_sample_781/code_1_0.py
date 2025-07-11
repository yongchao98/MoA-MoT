import math

# The problem asks for the largest number n of subcontinua in a decomposition
# of a specific type of topological space X.
# The space X has k=5 special points, and any subcontinuum containing m=3 of them must be X itself.
# We are looking for the largest n such that X = A_1 U ... U A_n, where each A_i is a proper
# subcontinuum and has a part not contained in the union of the others.
#
# A theorem in continuum theory states that the maximum size n for such a decomposition
# is given by the binomial coefficient C(k, m-1).

# Number of special points
k = 5

# Number of points that force a subcontinuum to be the whole space
m = 3

# The formula is n = C(k, m-1)
n = math.comb(k, m - 1)

# The final equation is C(5, 2) = 10.
# We will print each number in this final equation as requested.
print(k)
print(m - 1)
print(n)