import math

# Given parameters
r = 3
b = 9
false_positive_rate = 0.05

# We need to solve the equation: P(s) = 1 - (1 - s^r)^b for s.
# The formula to calculate s directly is: s = (1 - (1 - P(s))^(1/b))^(1/r)
s = (1 - (1 - false_positive_rate)**(1/b))**(1/r)

# Print the context and the equation to be solved
print("The goal is to find the similarity threshold 's' that satisfies the equation for the S-curve.")
print(f"The equation with r={r}, b={b}, and a false positive rate of {false_positive_rate} is:")
print(f"{false_positive_rate} = 1 - (1 - s^{r})^{b}")

# Print the final result rounded to three decimal places
print(f"\nThe calculated similarity threshold 's' should be: {s:.3f}")