# Define the parameters from the problem description
# c: Confidence of the discriminatory pattern in the subgraph G
c = 0.85
# s: Support of the discriminatory pattern in the subgraph G
s = 0.12
# c_prime: Confidence of the anti-correlated pattern in the 2-hop neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# Step 1: Calculate the strength of the evidence for the violation.
# This is the product of the confidence and support of the initial finding.
evidence_for_violation = c * s

# Step 2: Calculate the strength of the structurally-adjusted counter-evidence.
# The raw counter-evidence (c_prime) is modulated by the bottleneck coefficient (beta).
# A strong bottleneck (low beta) provides a structural explanation for the discrepancy,
# thus reducing the weight of the counter-evidence.
adjusted_counter_evidence = c_prime * beta

# Step 3: Calculate the final evidence ratio (E_r).
# This ratio compares the evidence for the violation against the counter-evidence.
evidence_ratio = evidence_for_violation / adjusted_counter_evidence

# Step 4: Print the result, showing the formula and the final value.
# The print statement includes all the original numbers to show the calculation.
print(f"The minimum ratio of evidence (E_r) is calculated as (c * s) / (c' * β).")
print(f"E_r = ({c} * {s}) / ({c_prime} * {beta})")
print(f"E_r = {evidence_for_violation} / {adjusted_counter_evidence}")
print(f"The final calculated evidence ratio is: {evidence_ratio}")
