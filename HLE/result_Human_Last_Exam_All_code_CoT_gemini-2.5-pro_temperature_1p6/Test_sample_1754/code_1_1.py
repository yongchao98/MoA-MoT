import math

# Define the constants lambda_1 and lambda_2
sqrt_17 = math.sqrt(17)
lambda_1 = (1 + sqrt_17) / 2
lambda_2 = (1 - sqrt_17) / 2

x2_0 = 0

# Calculate the value of the expression
value = (
    (2/3 * lambda_1 * math.exp(lambda_1 / 2) - 1/3 * math.exp(lambda_1 / 2)) * x2_0
    - (2/3 * lambda_2 * math.exp(lambda_2 / 2))
    - (10/3 * math.exp(lambda_1 / 2))
)

# Print the equation with substituted values.
# Each number is a component of the final sum.
# Let's calculate each component separately for clarity.
c1 = 2/3 * lambda_1 * math.exp(lambda_1/2)
c2 = -1/3 * math.exp(lambda_1/2)
c3 = -(2/3) * lambda_2 * math.exp(lambda_2/2)
c4 = -(10/3) * math.exp(lambda_1/2)

print("The full expression to calculate is:")
print(f"((2/3)*{lambda_1:.4f}*exp({lambda_1/2:.4f}) - (1/3)*exp({lambda_1/2:.4f})) * {x2_0} - (2/3)*{lambda_2:.4f}*exp({lambda_2/2:.4f}) - (10/3)*exp({lambda_1/2:.4f})")
print("\nSubstituting the values, the equation becomes:")
print(f"({c1:.4f} + {c2:.4f}) * {x2_0} + ({c3:.4f}) + ({c4:.4f}) = {value:.4f}")
print("\nSince x2(0) = 0, the equation simplifies to:")
print(f"0 + {c3:.4f} + {c4:.4f} = {value:.4f}")
print(f"\nFinal Answer: {value}")