import math

# Define the constants from the problem
sqrt5 = math.sqrt(5)
lambda1 = (1 + sqrt5) / 2
lambda2 = (1 - sqrt5) / 2

# Define the components of the final expression
c1 = 1
c5 = 5
base = 10
power = 10

# Calculate the two main terms in the sum
term1 = base**(power * lambda1)
term2 = base**(power * lambda2)

# Calculate the final result
result = (term1 + term2) / sqrt5

# Output the equation with the computed values of the constants
print(f"The expression to calculate is derived to be:")
print(f"({c1}/sqrt({c5})) * ({base}**({power}*{lambda1}) + {base}**({power}*{lambda2}))")
print("\nWhich simplifies to:")
print(f"(1 / {sqrt5}) * ({term1} + {term2})")
print("\nThe final calculated value is:")
print(result)
