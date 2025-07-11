import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of the student's code ---

# Define the variables provided in the problem description
c = 0.85  # Confidence of the discriminatory pattern P in subgraph G
s = 0.12  # Support for the pattern P in subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# Step 1: Calculate the evidence for the violation.
# This is derived from the pattern found in the local subgraph G.
# We combine confidence and support to measure its strength and prevalence.
evidence_for_violation = c * s

# Step 2: Calculate the effective counter-evidence.
# This comes from the anti-correlated pattern in the 2-hop neighborhood.
# Its influence is discounted by the bottleneck coefficient (beta), which
# represents the structural isolation of the subgraph G.
discounted_counter_evidence = c_prime * beta

# Step 3: Calculate the minimum evidence ratio (E_r).
# This is the ratio of the evidence for the violation to the discounted counter-evidence.
# This ratio tells us how much stronger the primary evidence is compared to the
# structurally weakened contradictory evidence.
E_r = evidence_for_violation / discounted_counter_evidence

# Output the explanation and the final equation with all values.
print("To determine the minimum ratio of evidence (E_r), we compare the evidence for the violation with the discounted counter-evidence.")
print("\nThe strength of the evidence for the violation in the subgraph is calculated as:")
print(f"Evidence_For = c * s")
print(f"             = {c} * {s} = {evidence_for_violation}")

print("\nThe strength of the counter-evidence from the neighborhood is discounted by the bottleneck coefficient β:")
print(f"Discounted_Counter_Evidence = c' * β")
print(f"                            = {c_prime} * {beta} = {discounted_counter_evidence}")

print("\nThe final evidence ratio (E_r) is the ratio of these two values:")
print(f"E_r = Evidence_For / Discounted_Counter_Evidence")
print(f"E_r = {evidence_for_violation} / {discounted_counter_evidence}")
print(f"E_r = {E_r}")

# --- End of the student's code ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the final output to the user
print(output)

# Automatically extract the final numerical answer for the platform
final_answer = None
try:
    # Find the last line that starts with "E_r ="
    lines = output.strip().split('\n')
    for line in reversed(lines):
        if line.startswith("E_r ="):
            final_answer = float(line.split('=')[1].strip())
            break
except (IndexError, ValueError) as e:
    # Fallback in case of parsing error
    final_answer = E_r
finally:
    if final_answer is not None:
        print(f"\n<<<{final_answer:.3f}>>>")
