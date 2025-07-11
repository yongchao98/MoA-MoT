# Define the given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern in subgraph G
s = 0.12  # Support of the discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlation in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# The goal is to calculate E_r, the ratio of evidence confirming a true violation
# versus it being a statistical artifact.

# Step 1: Calculate the evidence for a true violation.
# The raw evidence strength from the local observation is c * s.
# The hourglass structure (bottleneck) discounts the reliability of local evidence.
# A lower beta means a stronger bottleneck, so we multiply by beta to get the
# structurally-discounted evidence for the violation.
evidence_for_violation = (c * s) * beta

# Step 2: Calculate the evidence for a statistical artifact.
# The evidence for an artifact is supported by two things:
# 1. The direct counter-evidence from the wider neighborhood (c_prime).
# 2. The structural property of the graph that makes artifacts likely. This is
#    quantified by the strength of the bottleneck, (1 - beta).
# We combine these to find the total evidence for the artifact hypothesis.
evidence_for_artifact = c_prime * (1 - beta)

# Step 3: Compute the final ratio of evidence, E_r.
# This ratio compares the strength of the two competing hypotheses.
E_r = evidence_for_violation / evidence_for_artifact

# Step 4: Output the results, showing the formula and the values used.
print("The problem is to find the ratio of evidence (E_r) for a true fairness violation versus a statistical artifact.")
print("The formula is derived by comparing the discounted evidence for the violation against the evidence for the artifact explanation.")
print("\n--- Calculation Steps ---")
print(f"1. Evidence for Violation = (c * s) * beta = ({c} * {s}) * {beta} = {evidence_for_violation:.5f}")
print(f"2. Evidence for Artifact = c' * (1 - beta) = {c_prime} * (1 - {beta}) = {evidence_for_artifact:.4f}")
print(f"3. Ratio of Evidence (E_r) = Evidence for Violation / Evidence for Artifact")

print("\n--- Final Equation ---")
print(f"E_r = (c * s * beta) / (c' * (1 - beta))")
print(f"E_r = ({c} * {s} * {beta}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {evidence_for_violation:.5f} / {evidence_for_artifact:.4f}")
print(f"\nThe minimum ratio of evidence (E_r) is: {E_r}")