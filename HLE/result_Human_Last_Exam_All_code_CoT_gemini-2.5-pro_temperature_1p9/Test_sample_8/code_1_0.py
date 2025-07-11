import sympy

# Define the symbol for alpha
alpha = sympy.Symbol('α')

# Initial conormal order of the function f
s = 2 + alpha

# Order of the elliptic differential operator P = -Box_g - sigma^2
m = 2

# According to the theorem of elliptic regularity for conormal distributions,
# the resulting function u = R(sigma)f = P^{-1}f has a conormal order of s + m.
final_order = s + m

# The resulting function space
result_space = f"A^({sympy.printing.latex(final_order)})(X)"

# Print the step-by-step reasoning
print("The problem is to determine the conormal space of R(σ)f.")
print("Let u = R(σ)f. This is the solution to the equation (-□_g - σ^2)u = f.")
print("\nStep 1: Identify the order of the operator.")
print(f"The operator P = -□_g - σ^2 is a second-order elliptic differential operator, so its order is m = {m}.")

print("\nStep 2: Identify the conormal order of the input function.")
print(f"The function f belongs to A^(2+α)(X), so its conormal order is s = {s}.")

print("\nStep 3: Apply the elliptic regularity theorem.")
print("The theorem states that the conormal order of u is the sum of the order of f and the order of the operator P.")
print(f"Resulting order = s + m")
print(f"                 = ({s}) + {m}")
print(f"                 = {final_order}")

print("\nConclusion:")
print(f"Therefore, the function R(σ)f belongs to the conormal space: {result_space}")
