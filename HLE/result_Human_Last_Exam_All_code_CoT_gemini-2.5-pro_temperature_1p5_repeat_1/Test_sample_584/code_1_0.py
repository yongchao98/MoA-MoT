import math

# Define the LSH parameters and the target false positive rate
r = 3
b = 9
fp_rate = 0.05

# The equation for the threshold 's' is derived from P(s) = fp_rate,
# where P(s) = 1 - (1 - s^r)^b.
# Solving for s, we get: s = (1 - (1 - fp_rate)^(1/b))^(1/r)

# Calculate the threshold value
threshold_value = (1 - (1 - fp_rate)**(1/b))**(1/r)

# Round the result to three decimal points
rounded_threshold = round(threshold_value, 3)

# Print the final equation with all the numbers
print("To find the similarity threshold 's' where the false positive rate is < 0.05, we solve the following equation:")
print(f"s = (1 - (1 - {fp_rate})**(1/{b}))**(1/{r})")

# Print the final calculated threshold
print(f"\nThe calculated similarity threshold is: {rounded_threshold}")