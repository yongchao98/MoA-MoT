import sys

# Given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern in the subgraph G
s = 0.12  # Support of the discriminatory pattern in the subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
beta = 0.23  # Hourglass structural principle bottleneck coefficient

# Step 1: Calculate the strength of the direct evidence for the fairness violation.
# This is quantified by combining the local confidence and support of the pattern.
evidence_for_violation = c * s

# Step 2: Calculate the strength of the confounding evidence from the neighborhood.
# This evidence is discounted by the bottleneck coefficient (beta). A low beta indicates
# a strong bottleneck, making neighborhood evidence less impactful on the subgraph.
confounding_evidence = c_prime * beta

# Step 3: Calculate the minimum evidence ratio (E_r).
# This ratio compares the direct evidence to the discounted confounding evidence.
# A value greater than 1.0 would suggest the evidence for a violation is stronger.
if confounding_evidence == 0:
    print("Cannot calculate ratio: Confounding evidence is zero.", file=sys.stderr)
    E_r = float('inf')
else:
    E_r = evidence_for_violation / confounding_evidence

# Step 4: Output the results as requested.
# The final output should include the equation with the numbers substituted in.
print("To determine the minimum evidence ratio (E_r), we use the following model:")
print("E_r = (Direct Evidence Strength) / (Discounted Confounding Evidence Strength)")
print("E_r = (c * s) / (c' * β)")
print("\nPlugging in the given values:")
print(f"E_r = ({c} * {s}) / ({c_prime} * {beta})")
print(f"E_r = {evidence_for_violation:.4f} / {confounding_evidence:.4f}")
print(f"\nThe calculated minimum ratio of evidence required (E_r) is: {E_r}")
