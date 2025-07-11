# Given parameters from the automated knowledge graph FAIRness auditing system
# Confidence of the discriminatory pattern P in subgraph G
c = 0.85
# Support of the discriminatory pattern P in subgraph G
s = 0.12
# Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# Bottleneck coefficient of the knowledge graph
beta = 0.23

# --- Plan ---
# 1. The evidence for the violation is based on the local pattern's confidence (c) and support (s).
#    A combined score for this local evidence is c * s.
# 2. The hourglass structure with a bottleneck (beta) makes local patterns less reliable.
#    Therefore, the local evidence score is discounted by the bottleneck coefficient.
#    Adjusted evidence for violation = c * s * beta
# 3. The evidence for the pattern being a statistical artifact is the confidence of the
#    counter-evidence found in the wider neighborhood (c_prime).
# 4. The final evidence ratio (E_r) is the ratio of the adjusted evidence for the
#    violation to the evidence for the artifact.
#    Formula: E_r = (c * s * beta) / c_prime

# Calculate the numerator: the adjusted evidence for the violation
adjusted_violation_evidence = c * s * beta

# Calculate the final evidence ratio (E_r)
E_r = adjusted_violation_evidence / c_prime

# As requested, output all the numbers in the final equation
# The result shows a very low ratio, suggesting the evidence for the pattern being an artifact
# is much stronger than the evidence for it being a true violation.
print(f"Calculation: ({c} * {s} * {beta}) / {c_prime}")
print(f"Evidence Ratio (E_r): {E_r}")