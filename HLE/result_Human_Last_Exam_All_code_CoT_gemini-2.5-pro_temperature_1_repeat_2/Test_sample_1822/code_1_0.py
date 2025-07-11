import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of the logic ---

# Define the given parameters from the problem description
c = 0.85  # Confidence of the discriminatory pattern P in subgraph G
c_prime = 0.78  # Confidence of the anti-correlation in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# The plan is to calculate the evidence ratio (E_r) which is the ratio of
# the primary evidence strength to the confounding evidence strength.

# 1. Primary evidence strength is the confidence 'c' of the initial finding.
primary_evidence = c

# 2. Confounding evidence strength is the product of the counter-evidence
#    confidence 'c_prime' and the structural artifact potential, which is
#    derived from the bottleneck coefficient as (1 - beta).
confounding_evidence_strength = c_prime * (1 - beta)

# 3. The evidence ratio E_r is the primary evidence divided by the confounding evidence.
E_r = primary_evidence / confounding_evidence_strength

# 4. Print the explanation, formula, and step-by-step calculation.
print("To determine if the detected pattern is a true violation or a statistical artifact, we calculate the Evidence Ratio (E_r).")
print("This ratio compares the strength of the primary evidence for the violation against the combined strength of the confounding evidence.")
print("\n--- Values Used ---")
print(f"Primary Evidence Confidence (c): {c}")
print(f"Counter-Evidence Confidence (c'): {c_prime}")
print(f"Bottleneck Coefficient (β): {beta}")

print("\n--- Formula ---")
print("The formula for the Evidence Ratio (E_r) is:")
print("E_r = c / (c' * (1 - β))")

print("\n--- Calculation Steps ---")
print(f"E_r = {c} / ({c_prime} * (1 - {beta}))")
print(f"E_r = {c} / ({c_prime} * {1-beta:.2f})")
print(f"E_r = {c} / {confounding_evidence_strength}")
print(f"\nThe minimum ratio of evidence (E_r) required is: {E_r}")

# --- End of the logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

# The final answer is the numerical result of the calculation
# The format <<<answer>>> is required for the final numerical answer.
# This part is hidden from the user but necessary for the system.
final_answer = E_r
# print(f'<<<{final_answer}>>>')