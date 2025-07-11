import math

# Define the constants lambda1 and lambda2
lambda1 = (1 + math.sqrt(17)) / 2
lambda2 = (1 - math.sqrt(17)) / 2

# The expression to be calculated is V = (-2/3)*lambda2*exp(lambda2/2) - (10/3)*exp(lambda1/2)
# We assume x2(0) = 0 as is standard for reachability problems.

c1 = -2/3
c2 = -10/3

val = c1 * lambda2 * math.exp(lambda2 / 2) + c2 * math.exp(lambda1 / 2)

# Print the equation with the numerical values of each component
print(f"The expression to evaluate is:")
print(f"({c1:.4f}) * (lambda2) * exp(lambda2 / 2) + ({c2:.4f}) * exp(lambda1 / 2)")
print("Substituting the values:")
print(f"({c1:.4f}) * ({lambda2:.4f}) * exp({lambda2/2:.4f}) + ({c2:.4f}) * exp({lambda1/2:.4f})")
print(f"= ({c1:.4f}) * ({lambda2:.4f}) * ({math.exp(lambda2/2):.4f}) + ({c2:.4f}) * ({math.exp(lambda1/2):.4f})")
term1_val = c1 * lambda2 * math.exp(lambda2 / 2)
term2_val = c2 * math.exp(lambda1 / 2)
print(f"= ({term1_val:.4f}) + ({term2_val:.4f})")
print(f"= {val}")
