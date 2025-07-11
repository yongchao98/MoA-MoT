# Define the given parameters from the problem description
# c: Confidence of the initial discriminatory pattern
c = 0.85
# s: Support of the initial discriminatory pattern
s = 0.12
# c_prime: Confidence of the conflicting anti-correlation pattern in the neighborhood
c_prime = 0.78
# beta: Bottleneck coefficient of the knowledge graph
beta = 0.23

# --- Step 1: Calculate the evidence for the violation ---
# This is modeled as the product of the confidence and support of the local pattern.
# It represents the total weight of the observed discriminatory evidence.
evidence_for_violation = c * s

# --- Step 2: Calculate the evidence for it being a statistical artifact ---
# This is modeled by the confidence of the conflicting pattern, weighted by a
# structural factor derived from the bottleneck coefficient.
# (1 - beta) represents the structural propensity for artifacts: a stronger
# bottleneck (lower beta) results in a higher propensity.
evidence_for_artifact = c_prime * (1 - beta)

# --- Step 3: Calculate the final evidence ratio (E_r) ---
# E_r is the ratio of the pro-violation evidence to the pro-artifact evidence.
E_r = evidence_for_violation / evidence_for_artifact

# --- Step 4: Output the results ---
# The final output will present the equation with the plugged-in values and the final result.
print(f"The evidence ratio (E_r) is calculated using the formula:")
print(f"E_r = (c * s) / (c' * (1 - β))")
print(f"Plugging in the values:")
print(f"E_r = ({c} * {s}) / ({c_prime} * (1 - {beta}))")
print(f"E_r = {evidence_for_violation} / ({c_prime} * {1-beta})")
print(f"E_r = {evidence_for_violation} / {evidence_for_artifact}")
print(f"Final Result: E_r = {E_r}")

# Final Answer
print(f"\n<<<The minimum ratio of evidence required is {E_r}>>>")