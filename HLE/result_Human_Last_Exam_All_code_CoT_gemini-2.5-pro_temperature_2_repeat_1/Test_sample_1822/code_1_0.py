import sys

# Given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern in subgraph G
s = 0.12  # Support of the discriminatory pattern in subgraph G
c_prime = 0.78  # Confidence of the anti-correlation in the 2-hop neighborhood
beta = 0.23   # Bottleneck coefficient of the knowledge graph

# Step 1: Quantify the evidence for the violation based on local observations.
# This is a combination of the pattern's confidence and support.
evidence_violation = c * s

# Step 2: Quantify the counter-evidence from the neighborhood, adjusted for graph structure.
# The bottleneck coefficient (beta) modulates the influence of the neighborhood's pattern.
# A small beta means the neighborhood is less influential/relevant to the local subgraph.
evidence_artifact = c_prime * beta

# Step 3: Calculate the ratio of evidence (E_r).
# This ratio compares the strength of the evidence for a true violation against
# the strength of the evidence for it being a statistical artifact.
if evidence_artifact == 0:
    # Handle division by zero case, though not expected with given values.
    print("Error: Evidence for a statistical artifact is zero, ratio cannot be computed.")
    sys.exit()

E_r = evidence_violation / evidence_artifact

# Step 4: Output the results, explaining each component of the calculation.
print("This script calculates the evidence ratio (E_r) to assess a potential fairness violation.")
print("------------------------------------------------------------------------------------")
print(f"Given values:")
print(f"  - Local pattern confidence (c) = {c}")
print(f"  - Local pattern support (s) = {s}")
print(f"  - Neighborhood anti-correlation confidence (c') = {c_prime}")
print(f"  - Bottleneck coefficient (β) = {beta}\n")

print("Calculation Steps:")
print(f"1. The evidence for a fairness violation is calculated from local data: c * s")
print(f"   Evidence = {c} * {s} = {evidence_violation:.4f}\n")

print("2. The counter-evidence (for a statistical artifact) is taken from the neighborhood and adjusted by the graph's bottleneck structure: c' * β")
print(f"   Counter-Evidence = {c_prime} * {beta} = {evidence_artifact:.4f}\n")

print("3. The final ratio of evidence E_r is the evidence for violation divided by the counter-evidence:")
print(f"   E_r = (c * s) / (c' * β)")
print(f"   E_r = ({c} * {s}) / ({c_prime} * {beta})")
print(f"   E_r = {evidence_violation:.4f} / {evidence_artifact:.4f}")
print(f"   E_r = {E_r:.4f}\n")

print(f"The minimum ratio of evidence (E_r) is {E_r:.4f}.")

print(f"<<<{E_r:.4f}>>>")