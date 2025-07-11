import math

# Given parameters for the S-curve
r = 3
b = 9
# The target false positive rate (probability)
p = 0.05

# The formula to find the similarity threshold 's' is derived by solving p = 1 - (1 - s^r)^b for s.
# The resulting formula is: s = (1 - (1 - p)^(1/b))^(1/r)

# Calculate the threshold s
s_threshold = (1 - (1 - p)**(1/b))**(1/r)

# Print the equation with the given values substituted
print(f"The equation to solve for the threshold 's' is:")
print(f"s = (1 - (1 - {p})^(1/{b}))^(1/{r})")

# Print the final calculated threshold rounded to three decimal places
print(f"\nThe similarity threshold 's' is: {s_threshold:.3f}")