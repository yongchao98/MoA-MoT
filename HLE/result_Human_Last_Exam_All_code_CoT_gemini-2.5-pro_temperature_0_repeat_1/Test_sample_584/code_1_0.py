import math

# Given parameters for the S-curve
r = 3
b = 9
# The target false positive rate
fp_rate = 0.05

# The S-curve probability is given by P(s) = 1 - (1 - s^r)^b.
# We need to find the threshold 's' for which P(s) < fp_rate.
# The inequality is: 1 - (1 - s^r)^b < fp_rate
#
# Solving for s:
# (1 - s^r)^b > 1 - fp_rate
# 1 - s^r > (1 - fp_rate)^(1/b)
# s^r < 1 - (1 - fp_rate)^(1/b)
# s < (1 - (1 - fp_rate)^(1/b))^(1/r)

# Calculate the threshold value
threshold = (1 - (1 - fp_rate)**(1/b))**(1/r)

# Print the explanation and the final equation with numbers
print("To find the threshold 's' where the false positive rate is less than 0.05, we solve the following inequality:")
print(f"1 - (1 - s^{r})^{b} < {fp_rate}")
print("\nSolving for 's' gives us the expression for the threshold:")
print(f"s < (1 - (1 - {fp_rate})^(1/{b}))^(1/{r})")
print(f"\nThe calculated threshold is s < {threshold:.3f}")
print("\nTherefore, the threshold should be 0.178.")
