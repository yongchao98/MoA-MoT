import sympy as sp

# Define the symbolic variable s
s = sp.Symbol('s')

# A left coprime factorization of H(s) is H(s) = D(s)^-1 * N(s).
# The calculated matrices D(s) and N(s) are defined below.

# Denominator matrix D(s)
D = sp.Matrix([
    [s + 2, s - 1],
    [s + 1, 0]
])

# Numerator matrix N(s)
N = sp.Matrix([
    [s,     s + 2],
    [s - 1, s + 1]
])

# Print the resulting matrices
print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
print("\n")
print("D(s) =")
sp.pprint(D, use_unicode=True)
print("\n")
print("N(s) =")
sp.pprint(N, use_unicode=True)

# Verification of the factorization (optional)
# This part calculates D_inv * N to show it equals H(s).
print("\n" + "="*40)
print("Verification: Calculating D(s)^-1 * N(s)")
print("="*40)

# Calculate H(s) from D(s) and N(s)
H_calculated = sp.simplify(D.inv() * N)

print("H(s) = D(s)^-1 * N(s) =")
sp.pprint(H_calculated, use_unicode=True)
