import math

# The length of segment AB is given.
ab_length = 18468

# From the geometric analysis, we derived the formula for the length of CP.
# The relationship is CP = AB * (sqrt(5) + 1) / 2, where (sqrt(5) + 1) / 2 is the golden ratio.
sqrt_5 = math.sqrt(5)
golden_ratio = (sqrt_5 + 1) / 2

# Calculate the length of CP.
cp_length = ab_length * golden_ratio

# Print the breakdown of the calculation as requested.
print(f"Given the length of AB = {ab_length}")
print("The formula to calculate CP is: CP = AB * (sqrt(5) + 1) / 2")
print("\nCalculating the values:")
print(f"sqrt(5) = {sqrt_5}")
print(f"The equation is: CP = {ab_length} * (({sqrt_5} + 1) / 2)")
print(f"CP = {ab_length} * {golden_ratio}")
print(f"\nThe final length of CP is: {cp_length:.4f}")

# The final answer in the required format
# The value is rounded to 4 decimal points for the final output.
final_answer = round(cp_length, 4)
# print(f"<<<{final_answer}>>>")