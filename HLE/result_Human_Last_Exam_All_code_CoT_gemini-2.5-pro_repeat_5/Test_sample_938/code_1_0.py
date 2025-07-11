import math

# The problem simplifies to finding the limit L = ln(2 + sqrt(7)).
# We then need to compute the integer part of 10000 * L.

# Define the numbers in the final calculation
ten_thousand = 10000
two = 2
seven = 7

# Calculate the dominant root R
R = two + math.sqrt(seven)

# Calculate the limit L
L = math.log(R)

# Calculate the final value
final_value = ten_thousand * L

# Get the integer part of the final value
integer_part = int(final_value)

print(f"The dominant root R is 2 + sqrt(7) ≈ {R}")
print(f"The limit L is ln(R) ≈ {L}")
print(f"The expression to calculate is: {ten_thousand} * L")
print(f"The value of the expression is ≈ {final_value}")
print(f"The integer part of the value is: {integer_part}")
