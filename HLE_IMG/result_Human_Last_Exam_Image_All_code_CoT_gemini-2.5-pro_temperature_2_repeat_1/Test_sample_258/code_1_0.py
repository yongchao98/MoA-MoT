import math

# The length of the segment AB is given.
ab_length = 18468

# From the geometric analysis, the length of CP is found to be the length of AB
# multiplied by the golden ratio, phi.
# The formula is: CP = AB * (1 + sqrt(5)) / 2.

# Let's define the numbers used in the final equation.
one = 1
five = 5
two = 2

# Calculate the length of CP using the formula.
cp_length = ab_length * (one + math.sqrt(five)) / two

# As requested, we show the final equation with all the numbers
# and the result rounded to 4 decimal places.
print(f"Equation: CP = AB * (1 + sqrt(5)) / 2")
print(f"Calculation: CP = {ab_length} * ({one} + sqrt({five})) / {two}")
print(f"Result: CP = {cp_length:.4f}")
