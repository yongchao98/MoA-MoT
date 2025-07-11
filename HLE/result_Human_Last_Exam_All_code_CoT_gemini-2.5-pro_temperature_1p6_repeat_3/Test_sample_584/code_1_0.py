import math

# Given parameters for the S-curve
r = 3
b = 9
# The target false positive rate
false_positive_rate = 0.05

# The S-curve formula is P(s) = 1 - (1 - s^r)^b
# We want to find the similarity threshold 's' where P(s) < 0.05.
# To find the boundary, we solve the equation P(s) = 0.05 for s.
# 0.05 = 1 - (1 - s^r)^b
#
# Rearranging the formula to solve for s:
# s = (1 - (1 - false_positive_rate)^(1/b))^(1/r)

# Perform the calculation
s_threshold = (1 - (1 - false_positive_rate)**(1/b))**(1/r)

# The user wants each number in the final equation to be printed.
print("To find the threshold 's', we solve the following equation:")
print(f"s = (1 - (1 - {false_positive_rate})**(1/{b}))**(1/{r})")

print("\nSubstituting the numbers into the equation:")
val_inside_power = 1 - false_positive_rate
print(f"s = (1 - ({val_inside_power})**(1/{b}))**(1/{r})")

# Final Answer
rounded_s = round(s_threshold, 3)
print(f"\nThe calculated similarity threshold s is: {rounded_s}")

<<<0.179>>>