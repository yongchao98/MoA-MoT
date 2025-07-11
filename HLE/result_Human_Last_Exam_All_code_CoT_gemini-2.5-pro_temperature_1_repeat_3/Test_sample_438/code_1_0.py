import sys
import io

# Let's set a hypothetical baseline probability for the disease.
# This represents the risk with a normal, non-predisposing HLA variant.
# We'll assume a baseline risk of 0.05%
baseline_probability = 0.0005

# The problem states the variant increases the probability by a 1000-fold factor.
increase_factor = 1000

# Calculate the new, higher probability with the risk-associated HLA variant.
new_probability = baseline_probability * increase_factor

# --- Output Section ---
print("This script calculates the change in disease risk based on an HLA variant.")
print("-" * 60)

# Print each number that will be used in the final equation
print(f"Assumed baseline probability of disease: {baseline_probability}")
print(f"Fold increase due to the HLA variant: {increase_factor}")
print("-" * 60)

# Print the final equation showing all the numbers
print("The new risk probability is calculated as follows:")
print(f"New Probability = Baseline Probability * Increase Factor")
print(f"Resulting Equation: {new_probability} = {baseline_probability} * {increase_factor}")
print("-" * 60)

# Display the result in a user-friendly format
print(f"The new probability is {new_probability}, which is equivalent to a {new_probability:.0%} risk.")
print("This demonstrates a significant increase in disease risk.")
