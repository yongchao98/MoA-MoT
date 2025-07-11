from fractions import Fraction

# We are solving for x(1) using the equation:
# -1/x(1) = (1^2)/2 - 1 + 1/8

# Define the terms of the equation using the Fraction class for precision
term1 = Fraction(1**2, 2)
term2 = Fraction(1, 1)
term3 = Fraction(1, 8)

# The right hand side (RHS) of the equation
rhs = term1 - term2 + term3

# Solve for x(1)
# -1/x(1) = rhs  => x(1) = -1/rhs
x1 = -1 / rhs

print("To find x(1), we solve the equation: -1/x(1) = (1^2)/2 - 1 + 1/8")
print(f"The equation with its numerical values is: -1/x(1) = {term1} - {term2} + {term3}")
print(f"The right hand side evaluates to: {rhs}")
print(f"So, -1/x(1) = {rhs}")
print(f"This gives x(1) = -1 / ({rhs})")
print(f"The final value is x(1) = {x1}")
