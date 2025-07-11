import sympy

# Define the symbolic variable 's'
s = sympy.symbols('s')

# Define the left coprime matrices D(s) and N(s)
D = sympy.Matrix([
    [s + 1, 0],
    [1, s - 1]
])

N = sympy.Matrix([
    [s - 1, s + 1],
    [1, 1]
])

# Print the results in a clear format
print("The left coprime factorization is H(s) = D(s)^-1 * N(s)")
print("\nWhere D(s) is:")
sympy.pretty_print(D)

print("\nAnd N(s) is:")
sympy.pretty_print(N)

# Verify the result by calculating D^-1 * N
H_calculated = D.inv() * N
H_simplified = sympy.simplify(H_calculated)

print("\nVerification: D(s)^-1 * N(s) simplifies to:")
sympy.pretty_print(H_simplified)

# Also printing the original H(s) for comparison
H_original = sympy.Matrix([
    [(s - 1)/(s + 1), 1],
    [2/(s**2 - 1), 0]
])
print("\nOriginal H(s):")
sympy.pretty_print(H_original)
