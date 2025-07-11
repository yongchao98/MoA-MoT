import math

# Step 1: Define constants and target time
t_val = 2 * math.sqrt(3)
t_squared = t_val**2
print(f"The time is t = 2 * sqrt(3), so t^2 = {t_squared}")

# Step 2: The trajectory equation simplifies to x(t) = x(0) - 3.
# The value we need to compute is x(t) = (18 + 6*sqrt(3))^(1/3) + (18 - 6*sqrt(3))^(1/3)
print("\nThe final position is given by the expression: (18 + 6*sqrt(3))^(1/3) + (18 - 6*sqrt(3))^(1/3)")

# Step 3: Calculate the terms inside the cube roots
sqrt_3 = math.sqrt(3)
term1_inner = 18 + 6 * sqrt_3
term2_inner = 18 - 6 * sqrt_3

print(f"\nLet a = 18 + 6*sqrt(3). Calculated value: a = {term1_inner}")
print(f"Let b = 18 - 6*sqrt(3). Calculated value: b = {term2_inner}")

# Step 4: Calculate the cube roots
term1_cbrt = term1_inner**(1/3)
term2_cbrt = term2_inner**(1/3)

print(f"\nLet c = a^(1/3). Calculated value: c = {term1_cbrt}")
print(f"Let d = b^(1/3). Calculated value: d = {term2_cbrt}")

# Step 5: Sum the terms to get the final answer
final_position = term1_cbrt + term2_cbrt

print(f"\nThe final position is c + d.")
print(f"Final equation: {term1_cbrt} + {term2_cbrt} = {final_position}")