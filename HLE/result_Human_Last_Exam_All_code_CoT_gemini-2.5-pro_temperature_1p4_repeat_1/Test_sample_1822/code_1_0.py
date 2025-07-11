# Define the given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern in the subgraph
s = 0.12  # Support of the discriminatory pattern in the subgraph
c_prime = 0.78  # Confidence of the anti-correlated pattern in the neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# The minimum ratio of evidence (E_r) is calculated by comparing the strength of the
# initial finding against the strength of the counter-evidence, adjusted for the
# structural properties of the graph.
#
# Formula: E_r = (c * s) / (c' * (1 - beta))
# - Numerator (c * s): Represents the strength of the initial evidence for the discriminatory pattern.
# - Denominator (c' * (1 - beta)): Represents the effective strength of the counter-evidence.
#   The (1 - beta) term adjusts the counter-evidence confidence based on the graph's
#   structural bottleneck, which can isolate the subgraph from its neighborhood.

# Calculate the numerator (initial evidence strength)
initial_evidence_strength = c * s

# Calculate the denominator (effective counter-evidence strength)
effective_counter_evidence_strength = c_prime * (1 - beta)

# Calculate the final minimum ratio of evidence (E_r)
E_r = initial_evidence_strength / effective_counter_evidence_strength

# Print the step-by-step calculation and the final equation
print("--- FAIRness Violation Evidence Calculation ---")
print(f"Initial confidence (c) = {c}")
print(f"Initial support (s) = {s}")
print(f"Counter-evidence confidence (c') = {c_prime}")
print(f"Bottleneck coefficient (β) = {beta}\n")

print("Calculating the Minimum Evidence Ratio (E_r)...")
print("Formula: E_r = (c * s) / (c' * (1 - β))\n")

print("Final Equation with values:")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {initial_evidence_strength} / {effective_counter_evidence_strength}\n")

print(f"Result: The minimum ratio of evidence required (E_r) is {E_r}")

print(f"\n<<<{E_r}>>>")