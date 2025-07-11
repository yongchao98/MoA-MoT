import sympy

# Define the symbol and the function
x = sympy.Symbol('x')
y1 = 2 * x**-3

# Define the integration bounds
a = 4
b = 8

# Calculate the definite integral
# Indefinite integral first
integral_y1 = sympy.integrate(y1, x)
# Evaluate at the bounds
val_at_b = integral_y1.subs(x, b)
val_at_a = integral_y1.subs(x, a)
result = val_at_b - val_at_a

# Output the equation and the result
print(f"The integral is calculated as: -(1/{b}**2) - (-(1/{a}**2))")
print(f"= -1/{b**2} + 1/{a**2}")
print(f"= -1/64 + 4/64")
print(f"= {result}")
