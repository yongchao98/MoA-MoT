import math

# Given parameters for the S-curve
r = 3
b = 9
# The desired false positive rate threshold
fp_rate_threshold = 0.05

# The S-curve equation is P(s) = 1 - (1 - s^r)^b
# We need to find the similarity threshold 's' by solving the equation
# for when the probability equals the false positive rate.
# Equation: 1 - (1 - s^r)^b = fp_rate_threshold

# To solve for s, we rearrange the formula:
# s = (1 - (1 - fp_rate_threshold)^(1/b))^(1/r)

# Perform the calculation with the given values
s = (1 - (1 - fp_rate_threshold)**(1/b))**(1/r)

# Print the context and the equation being solved
print("We are solving for the similarity threshold 's' in the S-curve equation P(s) = 1 - (1 - s^r)^b.")
print("The goal is to find the value of 's' where the false positive rate is 0.05.")
print("\nGiven parameters:")
print(f"r = {r}")
print(f"b = {b}")
print(f"False Positive Rate = {fp_rate_threshold}")

print("\nThe equation to solve is:")
# Outputting each number in the final equation as requested
print(f"1 - (1 - s^{r})^{b} = {fp_rate_threshold}")


# Round the final result to three decimal points
s_rounded = round(s, 3)

print(f"\nThe calculated similarity threshold 's' is: {s_rounded}")
