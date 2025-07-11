import math

# Given parameters for the S-curve
r = 3
b = 9
# The target false positive rate
fp_rate = 0.05

# The equation to solve for the similarity threshold 's' is derived from:
# fp_rate = 1 - (1 - s^r)^b
#
# Rearranging the equation to solve for s:
# s = (1 - (1 - fp_rate)**(1/b))**(1/r)

# Calculate the threshold s
s_threshold = (1 - (1 - fp_rate)**(1/b))**(1/r)

# As requested, printing the final equation with each number
print("To find the threshold 's', we solve the following equation:")
print(f"s = (1 - (1 - {fp_rate})**(1/{b}))**(1/{r})")

# Print the final answer rounded to three decimal points
print(f"\nThe similarity threshold 's' should be: {s_threshold:.3f}")