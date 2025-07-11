import sympy

# Define the symbolic variable 's'
s = sympy.Symbol('s')

# Define the matrices D(s) and N(s) based on our derivation
D = sympy.Matrix([
    [1, s - 1],
    [0, s**2 - 1]
])

N = sympy.Matrix([
    [1, 1],
    [2, 0]
])

# Print the results in a clear format
print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
print("\nD(s) =")
sympy.pretty_print(D)

print("\nN(s) =")
sympy.pretty_print(N)

# The final answer in the required format
# Let's represent the matrices D and N row by row
# D = [[d11, d12], [d21, d22]]
# N = [[n11, n12], [n21, n22]]
d11, d12 = D.row(0)
d21, d22 = D.row(1)
n11, n12 = N.row(0)
n21, n22 = N.row(1)

print(f"""
Final Answer:
D(s) = 
| {d11}    {d12} |
| {d21}    {d22} |

N(s) = 
| {n11}    {n12} |
| {n21}     {n22} |
""")

# The final output is requested in a specific format " <<<answer content>>> ".
# I will output the matrices in a list-of-lists format.
# D(s) = [[1, s-1], [0, s**2-1]]
# N(s) = [[1, 1], [2, 0]]
answer_D = [[str(d11), str(d12)], [str(d21), str(d22)]]
answer_N = [[str(n11), str(n12)], [str(n21), str(n22)]]
answer_string = f"D(s) = {answer_D}, N(s) = {answer_N}"

# It's difficult to represent this structure directly in the required format.
# I will print the components of the D and N matrices.
final_answer_str = f"D(s) = [[{d11}, {d12}], [{d21}, {d22}]], N(s) = [[{n11}, {n12}], [{n21}, {n22}]]"

# Since a single simple answer is preferred, let's format it nicely.
# Let's provide the final result D(s) and N(s) in a structured way.
print("Final Answer in <<<...>>> format:")
print("D(s) matrix, row 1: [1, s - 1]")
print("D(s) matrix, row 2: [0, s**2 - 1]")
print("N(s) matrix, row 1: [1, 1]")
print("N(s) matrix, row 2: [2, 0]")
<<<D(s) = [[1, s-1], [0, s**2-1]], N(s) = [[1, 1], [2, 0]]>>>