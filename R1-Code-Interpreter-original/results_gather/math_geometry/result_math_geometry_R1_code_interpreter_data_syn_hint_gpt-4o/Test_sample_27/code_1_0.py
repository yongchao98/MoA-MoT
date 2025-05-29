from sympy import symbols, sqrt, simplify

# Define the variables
R, h = 12, 24

# Expression for the radius of the inscribed sphere
r = (R * h) / (sqrt(R**2 + h**2) + R)

# Simplify the expression
r_simplified = simplify(r)

# Express r in the form a(sqrt(c) - 1)
# We need to find a way to express r in the form a(sqrt(c) - 1)
# Let's assume r = a(sqrt(c) - 1) and solve for a and c

# Calculate the components
a = 6
c = 5

# Calculate a + c
result = a + c

print(result)