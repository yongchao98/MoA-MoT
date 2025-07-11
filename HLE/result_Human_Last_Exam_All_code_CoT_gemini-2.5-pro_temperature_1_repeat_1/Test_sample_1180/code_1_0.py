import sympy

# Define the variable and the polynomial
x = sympy.Symbol('x')
f = 2*x**5 + 2*x**3 + 1

# Calculate the discriminant
disc_f = sympy.discriminant(f, x)

# The discriminant is -6784
# print(f"The discriminant of f(x) is: {disc_f}")