# Define the given parameters
c = 0.85  # Confidence of the discriminatory pattern in subgraph G
s = 0.12  # Support for the pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the neighborhood
beta = 0.23  # Bottleneck coefficient of the hourglass structure

# Step 1: Calculate the strength of the confirming evidence for the violation.
# This is modeled as the product of the local confidence and support.
confirming_evidence = c * s

# Step 2: Calculate the strength of the conflicting evidence.
# This is modeled as the product of the neighborhood confidence and a factor
# representing the likelihood of a structural artifact, (1 - beta).
# A strong bottleneck (low beta) increases this factor, making the local
# finding more suspect.
conflicting_evidence = c_prime * (1 - beta)

# Step 3: Calculate the minimum ratio of evidence (E_r).
# This ratio compares the confirming evidence to the conflicting evidence.
E_r = confirming_evidence / conflicting_evidence

# Output the explanation and the final equation with all values.
print("Calculating the minimum ratio of evidence (E_r)")
print("Formula: E_r = (c * s) / (c' * (1 - beta))")
print("\nPlugging in the values:")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {confirming_evidence} / ({c_prime} * {1-beta})")
print(f"E_r = {confirming_evidence} / {conflicting_evidence}")
print(f"\nFinal Evidence Ratio (E_r): {E_r}")
