import sympy

# Set up the symbolic variable
s = sympy.Symbol('s')

# Define the calculated left coprime matrices D(s) and N(s)
# based on the derivation above.
D = sympy.Matrix([
    [1, s - 1],
    [0, s**2 - 1]
])

N = sympy.Matrix([
    [1, 1],
    [2, 0]
])

# Define the original transfer function H(s) for verification
H = sympy.Matrix([
    [(s - 1) / (s + 1), 1],
    [2 / (s**2 - 1), 0]
])

# --- Output the results ---

print("A left coprime factorization for the transfer function H(s) is H(s) = D(s)^-1 * N(s).")
print("\nThe matrices D(s) and N(s) are found to be:")

print("\nD(s) =")
sympy.pretty_print(D)

print("\nN(s) =")
sympy.pretty_print(N)

# The final equation with each matrix element shown:
print("\nThe final equation is:")

# Pretty print the full equation H(s) = D(s)^-1 N(s)
# We manually construct the print output for better alignment and clarity.
h_str = sympy.printing.pretty(H, use_unicode=True).split('\n')
d_str = sympy.printing.pretty(D, use_unicode=True).split('\n')
n_str = sympy.printing.pretty(N, use_unicode=True).split('\n')

# Add the inverse symbol '⁻¹' to D(s)
mid_row = len(d_str) // 2
d_str[mid_row] += '⁻¹'

# Align the matrices for printing
max_height = max(len(h_str), len(d_str), len(n_str))
while len(h_str) < max_height:
    h_str.insert(len(h_str)//2, ' ' * len(h_str[0]))
while len(d_str) < max_height:
    d_str.insert(len(d_str)//2, ' ' * len(d_str[0]))
while len(n_str) < max_height:
    n_str.insert(len(n_str)//2, ' ' * len(n_str[0]))

for i in range(max_height):
    op1 = '  =  ' if i == max_height // 2 else '     '
    op2 = ' * ' if i == max_height // 2 else '   '
    print(h_str[i] + op1 + d_str[i] + op2 + n_str[i])

# Verification (optional, for confidence)
# D_inv_N = sympy.simplify(D.inv() * N)
# print(f"\nVerification successful: {D_inv_N == H}")
