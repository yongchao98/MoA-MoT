import math

# The problem is to find the maximum achievable ratio of bidirectional conical power
# to the intensity along a line. The final expression for this ratio is:
# R_max = 4 * pi * (4/3 - 7 / (6 * sqrt(2)))

# We will calculate the numerical value of this expression.
pi_val = math.pi
sqrt2_val = math.sqrt(2)

# Here are the values of the individual numbers and terms in the final equation.
term1 = 4 / 3
term2 = 7 / (6 * sqrt2_val)
factor = 4 * pi_val

# The final result is calculated by combining these terms.
max_ratio = factor * (term1 - term2)

print("The final equation for the ratio is: (4 * pi) * (4/3 - 7 / (6 * sqrt(2)))")
print(f"The value for the term 4/3 is: {term1}")
print(f"The value for the term 7 / (6 * sqrt(2)) is: {term2}")
print(f"The value for the factor 4 * pi is: {factor}")
print(f"The maximum achievable ratio is: {max_ratio}")