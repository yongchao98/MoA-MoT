import math

# Given eigenvalues
sqrt_17 = math.sqrt(17)
lambda1 = (1 + sqrt_17) / 2
lambda2 = (1 - sqrt_17) / 2

# Based on the problem analysis, we assume the initial condition x(0)=0, so x2(0)=0.
x2_0 = 0

# Calculate each term of the expression
# The expression is: (2/3*lambda1*e^(lambda1/2) - 1/3*e^(lambda1/2))*x2(0) - 2/3*lambda2*e^(lambda2/2) - 10/3*e^(lambda1/2)

term1_factor = (2/3 * lambda1 * math.exp(lambda1 / 2) - 1/3 * math.exp(lambda1 / 2))
term1 = term1_factor * x2_0

term2 = (2/3) * lambda2 * math.exp(lambda2 / 2)

term3 = (10/3) * math.exp(lambda1 / 2)

# Calculate the final result
result = term1 - term2 - term3

# Print the equation with the calculated values for each term as requested
print("The equation to be calculated is:")
print(f"({term1_factor}) * {x2_0} - ({term2}) - ({term3})")
print(f"= {term1} - ({term2}) - ({term3})")
print(f"= {result}")
