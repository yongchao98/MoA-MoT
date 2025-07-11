import math

# Given parameters
r = 3
b = 9
fp_rate = 0.05

# The S-curve equation is P(s) = 1 - (1 - s^r)^b
# We need to find the threshold s where P(s) = fp_rate
# So, we solve the equation: 1 - (1 - s^r)^b = fp_rate

# Rearranging the formula to solve for s:
# (1 - s^r)^b = 1 - fp_rate
# 1 - s^r = (1 - fp_rate)^(1/b)
# s^r = 1 - (1 - fp_rate)^(1/b)
# s = (1 - (1 - fp_rate)^(1/b))^(1/r)

s = (1 - (1 - fp_rate)**(1/b))**(1/r)

# Print the explanation and the final result
print(f"Given r = {r} and b = {b}, we want to find the similarity threshold 's' where the false positive rate is {fp_rate}.")
print(f"This requires solving the equation: 1 - (1 - s^{r})^{b} = {fp_rate}")
print(f"Plugging in the values, we get: 1 - (1 - s^{r})^{b} = {fp_rate}")
print(f"The threshold 's' is calculated as: (1 - (1 - {fp_rate})**(1/{b}))**(1/{r})")
print(f"The calculated threshold 's' is: {s:.3f}")

<<<0.178>>>