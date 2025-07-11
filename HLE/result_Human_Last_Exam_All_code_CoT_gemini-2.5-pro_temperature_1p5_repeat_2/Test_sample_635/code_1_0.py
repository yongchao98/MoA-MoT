import sympy

# Define the symbolic variable s
s = sympy.Symbol('s')

# Define the transfer function H(s)
H = sympy.Matrix([[(s - 1) / (s + 1), 1],
                  [2 / (s**2 - 1), 0]])

# Define the calculated left coprime factors D(s) and N(s)
D = sympy.Matrix([[1, s - 1],
                  [s + 1, 0]])

N = sympy.Matrix([[1, 1],
                  [s - 1, s + 1]])

# Verify the factorization H(s) = D(s)^-1 * N(s)
# Calculate D_inv_N = D(s)^-1 * N(s)
D_inv = D.inv()
D_inv_N = D_inv * N

# Simplify the resulting expression
D_inv_N_simplified = sympy.simplify(D_inv_N)

# Print the results
print("The transfer function H(s) is:")
sympy.pprint(H)
print("\n" + "="*30 + "\n")

print("The calculated left denominator matrix D(s) is:")
sympy.pprint(D)
print("\n" + "="*30 + "\n")

print("The calculated left numerator matrix N(s) is:")
sympy.pprint(N)
print("\n" + "="*30 + "\n")

print("Verification: D(s)^-1 * N(s) simplifies to:")
sympy.pprint(D_inv_N_simplified)
print("\n" + "="*30 + "\n")

# Check if the factorization is correct
is_correct = (sympy.simplify(H - D_inv_N) == sympy.zeros(2, 2))
print(f"Is the factorization D(s)^-1 * N(s) = H(s) correct? {is_correct}")
print("\nFinal Answer: The matrices D(s) and N(s) are given above.")

print("\nTo express the final answer in the required format, we print each element of D(s) and N(s).")

print("\nD(s) = ")
print("[[{}, {}],".format(D[0,0], D[0,1]))
print(" [{}, {}]]".format(D[1,0], D[1,1]))

print("\nN(s) = ")
print("[[{}, {}],".format(N[0,0], N[0,1]))
print(" [{}, {}]]".format(N[1,0], N[1,1]))
