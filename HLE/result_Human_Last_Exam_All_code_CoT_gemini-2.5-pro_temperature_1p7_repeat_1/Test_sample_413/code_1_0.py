import sympy

# Define the variable
x = sympy.Symbol('x')

# Define the coefficients based on the derivation
a = -sympy.Rational(3, 8)
b = -sympy.Rational(5, 8)
c = 0

# Define the function
f = x**3 + a*x**2 + b*x + c

# Calculate f(3)
val_f3 = f.subs(x, 3)

# Print the final equation
print(f"f(x) = x^3 + ({a})x^2 + ({b})x + {c}")
print(f"f(3) = 3^3 + ({a})*3^2 + ({b})*3 + {c}")
print(f"f(3) = 27 + ({a*9}) + ({b*3})")
print(f"f(3) = 27 - 27/8 - 15/8")
print(f"f(3) = 27 - 42/8")
print(f"f(3) = 27 - 21/4")
print(f"f(3) = 108/4 - 21/4")
print(f"f(3) = {val_f3}")