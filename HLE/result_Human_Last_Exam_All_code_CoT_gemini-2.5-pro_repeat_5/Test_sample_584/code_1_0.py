import math

# Given parameters for the S-curve 1-(1-s^r)^b
r = 3
b = 9

# The target false positive rate, which is the probability P(s)
fp_rate = 0.05

# We need to find the similarity threshold 's' that produces this rate.
# The equation to solve is: fp_rate = 1 - (1 - s^r)^b
#
# Rearranging the equation to solve for 's', we get:
# s = (1 - (1 - fp_rate)^(1/b))^(1/r)

# Now, we plug in the values and calculate 's'.
s_threshold = (1 - (1 - fp_rate)**(1/b))**(1/r)

# As requested, we first print the final equation with all numbers substituted.
print(f"To find the threshold 's' where the false positive rate is {fp_rate}, we solve the following equation:")
print(f"s = (1 - (1 - {fp_rate})**(1/{b}))**(1/{r})")

# Then, we print the calculated result rounded to three decimal points.
print(f"\nThe calculated similarity threshold should be: {s_threshold:.3f}")