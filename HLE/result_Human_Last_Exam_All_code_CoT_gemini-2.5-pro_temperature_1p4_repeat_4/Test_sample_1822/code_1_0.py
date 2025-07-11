# Define the given parameters
# Confidence of the discriminatory pattern in the subgraph G
c = 0.85
# Support of the discriminatory pattern in the subgraph G
s = 0.12
# Confidence of the anti-correlation in the 2-hop neighborhood
c_prime = 0.78
# Bottleneck coefficient of the knowledge graph
beta = 0.23

# Formulate the model for the evidence ratio (E_r)
# E_r = (Evidence for True Violation) / (Evidence for Statistical Artifact)

# Evidence for True Violation: The raw evidence (c * s) weighted by the reliability
# of local findings, which is proportional to the bottleneck coefficient (beta).
evidence_true = c * s * beta

# Evidence for Statistical Artifact: The counter-evidence (c_prime) weighted by the
# structural propensity for artifacts, modeled as (1 - beta).
evidence_artifact = c_prime * (1 - beta)

# Calculate the final evidence ratio
E_r = evidence_true / evidence_artifact

# Print the step-by-step calculation as requested
print("The formula for the minimum ratio of evidence (E_r) is:")
print("E_r = (c * s * β) / (c' * (1 - β))")
print("\nSubstituting the given values:")
print(f"E_r = ({c} * {s} * {beta}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = ({c * s * beta}) / ({c_prime * (1 - beta)})")
print(f"\nThe minimum ratio of evidence E_r is: {E_r}")
