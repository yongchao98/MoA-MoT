import math

# The equation for the generating amplitude c1, derived under the
# assumption of a corrected problem statement (standard van der Pol oscillator),
# is 2 * c1^2 = 4.

# Coefficients of the equation a * c1**p = b
a = 2
p = 2
b = 4

print("The equation for c1 is:")
# As requested, outputting each number in the final equation.
# The variable 'c1' is represented as a string.
print(a, "* c1**", p, "=", b)

# Solve for c1 > 0
c1_squared = b / a
first_positive_root_c1 = math.sqrt(c1_squared)

print("\nThe value of the first positive root c1 is:")
print(first_positive_root_c1)