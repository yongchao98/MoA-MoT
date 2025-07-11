import sympy

# Define the variable
t = sympy.symbols('t')

# Define the numerator and the terms for the denominator of the solution x(t)
numerator = -8
a = 4
b = -8
c = 1

# Value of t to evaluate at
t_val = 1

# Calculate the denominator at t=1
denominator = a * t_val**2 + b * t_val + c

# Calculate the result
result = sympy.Rational(numerator, denominator)

# Print the equation and the steps
print(f"From the derived solution x(t) = -8 / (4*t^2 - 8*t + 1)")
print(f"We substitute t = {t_val}:")
print(f"x({t_val}) = {numerator} / ({a}*({t_val})^2 - {abs(b)}*({t_val}) + {c})")
print(f"x({t_val}) = {numerator} / ({a*t_val**2} - {abs(b)*t_val} + {c})")
print(f"x({t_val}) = {numerator} / ({denominator})")
print(f"x({t_val}) = {result}")
