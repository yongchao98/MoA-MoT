import math

# Given parameters for the LSH S-curve
r = 3
b = 9

# The target false positive rate (FPR)
fpr = 0.05

# The S-curve equation gives the probability P(s) that a pair with similarity 's'
# becomes a candidate pair: P(s) = 1 - (1 - s^r)^b
# We want to find the threshold 's' for which P(s) equals our target FPR.

# Equation: fpr = 1 - (1 - s^r)^b
# Rearranging to solve for s:
# 1 - fpr = (1 - s^r)^b
# (1 - fpr)^(1/b) = 1 - s^r
# s^r = 1 - (1 - fpr)^(1/b)
# s = (1 - (1 - fpr)^(1/b))^(1/r)

# Calculate the threshold 's'
s_threshold = (1 - (1 - fpr)**(1/b))**(1/r)

print(f"Given r = {r}, b = {b}, and a target false positive rate of {fpr}.")
print("We need to find the similarity threshold 's' by solving the equation:")
print(f"{fpr} = 1 - (1 - s^{r})^{b}\n")
print("The final calculation is:")
print(f"s = (1 - (1 - {fpr})**(1/{b}))**(1/{r})")

print(f"\nThe threshold 's' is: {s_threshold:.3f}")