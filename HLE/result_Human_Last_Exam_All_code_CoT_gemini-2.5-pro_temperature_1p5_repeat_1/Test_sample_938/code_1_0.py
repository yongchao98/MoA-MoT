import math

# Step 1: Define the constants for the final calculation.
c1 = 10000
c2 = 2
c3 = 7

# Step 2: The limit is ln(alpha), where alpha is the dominant root of the characteristic equation.
# The dominant root is 2 + sqrt(7).
alpha = c2 + math.sqrt(c3)

# Step 3: Calculate the value of the limit expression.
limit_value = math.log(alpha)
final_value = c1 * limit_value

# Step 4: Output the result and the equation.
print(f"The limit is ln(2 + sqrt(7))")
print(f"The value to compute is {c1} * ln({c2} + sqrt({c3}))")
print(f"{c1} * {limit_value:.6f}... = {final_value:.6f}...")

# Step 5: Find and print the integer part.
integer_part = int(final_value)
print(f"The integer part of the value is {integer_part}.")

# Final answer in the specified format
# print(f"<<<{integer_part}>>>")