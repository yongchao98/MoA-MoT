import math

# Given parameters for the LSH S-curve
r = 3
b = 9
fp_rate = 0.05

# The S-curve equation is P(s) = 1 - (1 - s^r)^b
# We need to find the similarity 's' where the probability P(s) is equal to the desired false positive rate.
# So, we solve the equation: fp_rate = 1 - (1 - s^r)^b

# Algebraically solving for 's':
# 1 - fp_rate = (1 - s^r)^b
# (1 - fp_rate)^(1/b) = 1 - s^r
# s^r = 1 - (1 - fp_rate)^(1/b)
# s = (1 - (1 - fp_rate)**(1/b))**(1/r)

# Calculate the value of s
s = (1 - (1 - fp_rate)**(1/b))**(1/r)

# Print the final equation with the numbers substituted in
print("We are solving the following equation for 's':")
print(f"{fp_rate} = 1 - (1 - s^{r})^{b}")
print("\nSubstituting the given values r=3 and b=9:")
print(f"{fp_rate} = 1 - (1 - s^{{{r}}})^{{{b}}}")

# Print the final result rounded to three decimal places
print(f"\nThe similarity threshold 's' should be: {s:.3f}")