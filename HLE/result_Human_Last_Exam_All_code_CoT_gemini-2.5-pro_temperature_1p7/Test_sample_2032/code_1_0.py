import sympy

# Define the symbol for our variable x
x = sympy.Symbol('x')

# Define the derived functions for I(x) and J(x) for x in [0, 1/2].
# I(x) = E[Y | X_1=x]
# J(x) = E[Y^2 | X_1=x]
I_x = sympy.Rational(1, 2) - 4*x**3 + 8*x**4
J_x = sympy.Rational(3, 10) - 6*x**4 + sympy.Rational(64, 5)*x**5

# Calculate E[Y]
E_Y = sympy.Rational(1, 2)

# Calculate the definite integrals of I(x) and J(x) from 0 to 1/2
integral_I = sympy.integrate(I_x, (x, 0, sympy.Rational(1, 2)))
integral_J = sympy.integrate(J_x, (x, 0, sympy.Rational(1, 2)))

# Calculate E[Y^2] using the derived formula
E_Y_squared = sympy.Rational(1, 2) - 2 * integral_I + 2 * integral_J

# Calculate the Variance of Y
Var_Y = E_Y_squared - E_Y**2

# Print the results of each step
print(f"The expected value E[Y] = {E_Y}")
print(f"The integral of I(x) from 0 to 1/2 is: {integral_I}")
print(f"The integral of J(x) from 0 to 1/2 is: {integral_J}")
print(f"The second moment E[Y^2] = 1/2 - 2 * ({integral_I}) + 2 * ({integral_J}) = {E_Y_squared}")
print(f"The variance is calculated using the formula: Var(Y) = E[Y^2] - (E[Y])^2")
print(f"So, Var(Y) = {E_Y_squared} - ({E_Y})^2 = {E_Y_squared} - {E_Y**2} = {Var_Y}")
