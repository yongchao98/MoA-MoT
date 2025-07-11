# Define the parameters from the problem description
# Confidence of the initial discriminatory pattern
c = 0.85
# Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# Bottleneck coefficient based on the hourglass structural principle
beta = 0.23

# The minimum evidence ratio (E_r) is calculated by dividing the confidence
# of the initial pattern (evidence for) by the product of the counter-evidence
# confidence and the bottleneck coefficient (evidence against).
# A small bottleneck coefficient (beta) increases the denominator, reflecting
# that a higher standard of evidence is required to overcome structural uncertainty.

# Calculate the value for the denominator
denominator = c_prime * beta

# Calculate the final evidence ratio
E_r = c / denominator

# Print the final equation with all the numbers and the result
# The numbers in the equation are shown as required.
print("The minimum required ratio of evidence (E_r) is calculated as follows:")
print(f"E_r = c / (c' * β)")
print(f"E_r = {c} / ({c_prime} * {beta})")
print(f"E_r = {E_r:.3f}")
<<<4.738>>>