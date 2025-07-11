import sympy

# Set up the symbolic variable for 's'
s = sympy.Symbol('s')

# The final calculated left coprime matrices D(s) and N(s)
D_s = sympy.Matrix([
    [s + 1, 0],
    [1, s - 1]
])

N_s = sympy.Matrix([
    [s - 1, s + 1],
    [1, 1]
])

# Print the final result in the required format
print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by the polynomial matrices:")

print("\nD(s) =")
# Print each element of the D(s) matrix
print(f"[[{D_s[0, 0]}, {D_s[0, 1]}]")
print(f" [{D_s[1, 0]}, {D_s[1, 1]}]]")


print("\nN(s) =")
# Print each element of the N(s) matrix
print(f"[[{N_s[0, 0]}, {N_s[0, 1]}]")
print(f" [{N_s[1, 0]}, {N_s[1, 1]}]]")