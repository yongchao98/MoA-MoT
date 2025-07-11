import math

# Define the given parameters
r = 3
b = 9
false_positive_rate = 0.05

print(f"We are solving for the similarity threshold 's' using the S-curve formula: P = 1 - (1 - s^r)^b")
print(f"Given r = {r}, b = {b}, and a desired false positive rate (P) < {false_positive_rate}.")

# To find the threshold, we solve for 's' when P is exactly the false positive rate.
# The solved equation is: s = (1 - (1 - P)^(1/b))^(1/r)

# Now, we plug the numbers into the equation.
print("\nThe final equation with the given numbers is:")
print(f"s = (1 - (1 - {false_positive_rate})**(1/{b}))**(1/{r})")

# Perform the calculation
threshold = (1 - (1 - false_positive_rate)**(1/b))**(1/r)

# Round the result to three decimal points
rounded_threshold = round(threshold, 3)

print(f"\nThe calculated threshold where the false positive rate is less than 0.05 is: {rounded_threshold}")