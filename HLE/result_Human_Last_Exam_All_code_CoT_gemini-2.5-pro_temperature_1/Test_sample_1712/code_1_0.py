import math

# The value of the sum is given by the analytical expression:
# S = 1 - (pi^2)/12 + (ln(2))^2/2

# Calculate the value of each term in the expression.
term1 = 1.0
term2 = (math.pi**2) / 12
term3 = (math.log(2)**2) / 2

# Calculate the final result.
result = term1 - term2 + term3

# As requested, output the numbers in the final equation.
# The equation is: 1 - (pi^2)/12 + (ln(2))^2/2 = Result
print(f"{term1} - {term2} + {term3} = {result}")
