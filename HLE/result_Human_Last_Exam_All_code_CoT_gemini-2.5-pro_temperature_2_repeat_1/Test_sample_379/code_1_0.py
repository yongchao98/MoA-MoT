import math

# The length of the curve C is (4 * pi * r) / 3.
# The problem asks for the ratio of this length to (2 * pi * r).
# The value 'r' is a variable, but it will cancel out in the ratio.
# We can set r=1 for simplicity.
r = 1
length_C = (4 * math.pi * r) / 3
length_divisor = 2 * math.pi * r

# Calculate the ratio
ratio = length_C / length_divisor

# Convert the ratio to a percentage
percentage = ratio * 100

# Round to the nearest whole number
whole_number_percentage = round(percentage)

# Final Equation steps to be printed
# We need to output each number in the final equation
# The equation is: round(((4 * pi * r / 3) / (2 * pi * r)) * 100)
# Let's show the cancellation and calculation
# (4/3) / 2 * 100
a = 4
b = 3
c = 2
d = 100

final_val = (a/b) / c * d

print(f"The length of the curve C is (4 * pi * r) / 3.")
print(f"We divide this by 2 * pi * r.")
print(f"The ratio is ({a} * pi * r / {b}) / ({c} * pi * r).")
print(f"The terms 'pi' and 'r' cancel out, leaving ({a}/{b}) / {c}.")
print(f"This simplifies to {a}/({b}*{c}) = {a}/{b*c}, which is {round(a/(b*c), 3)}...")
print(f"To get a percentage, we multiply by {d}: {round(a/(b*c), 3)}... * {d} = {round(final_val, 1)}%")
print(f"Rounding to the nearest whole number, the answer is: {whole_number_percentage}")
