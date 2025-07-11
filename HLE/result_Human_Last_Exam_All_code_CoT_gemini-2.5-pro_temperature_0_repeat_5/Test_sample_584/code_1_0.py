import math

# Given parameters for the S-curve
r = 3
b = 9
# The target false positive rate
p = 0.05

# The S-curve equation is P(s) = 1 - (1 - s^r)^b
# We need to find the similarity threshold 's' where P(s) = p.
# To do this, we solve the equation for 's':
# s = (1 - (1 - p)^(1/b))^(1/r)

# Calculate the threshold
s_threshold = math.pow(1 - math.pow(1 - p, 1/b), 1/r)

# Print the equation with the given numbers
print("The equation to calculate the threshold 's' is:")
print(f"s = (1 - (1 - {p})^(1/{b}))^(1/{r})")

# Print the final calculated threshold, formatted to three decimal places
print(f"\nThe threshold 's' where the false positive rate is 0.05 is: {s_threshold:.3f}")
<<<0.179>>>