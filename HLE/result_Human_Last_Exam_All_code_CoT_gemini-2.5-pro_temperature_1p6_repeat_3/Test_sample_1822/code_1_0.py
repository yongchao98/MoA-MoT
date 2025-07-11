import sys

# Define the given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern in the subgraph G
s = 0.12  # Support of the discriminatory pattern in the subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# Step 1: Calculate the evidence FOR the fairness violation.
# This is modeled as the product of the confidence and support of the observed discriminatory pattern.
evidence_for_violation = c * s

# Step 2: Calculate the evidence AGAINST the violation (i.e., for it being a statistical artifact).
# This is modeled by combining the confidence of the counter-evidence (c') with the structural
# propensity for artifacts. A strong bottleneck (low beta) increases this propensity, which
# we model as (1 - beta).
evidence_against_violation = c_prime * (1 - beta)

# Step 3: Calculate the final evidence ratio (E_r).
# This ratio weighs the evidence supporting the violation against the counter-evidence.
# A value < 1 suggests the counter-evidence is stronger.
# A value > 1 suggests the evidence for the violation is stronger.
if evidence_against_violation == 0:
    # Handle division by zero case, although unlikely with the given inputs.
    print("Error: Counter-evidence calculation resulted in zero, cannot compute ratio.", file=sys.stderr)
    E_r = float('inf')
else:
    E_r = evidence_for_violation / evidence_against_violation

# Step 4: Print the equation with all the numbers and the final result.
print("The evidence ratio (E_r) is calculated by the formula:")
print(f"E_r = (c * s) / (c' * (1 - β))")
print("\nSubstituting the given values:")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {evidence_for_violation} / ({c_prime} * {1-beta})")
print(f"E_r = {evidence_for_violation} / {evidence_against_violation}")
print("\nFinal calculated ratio:")
print(f"E_r = {E_r}")

# The final answer in the requested format
# <<<answer>>>
final_answer = E_r
# <<<f"{final_answer}">>>