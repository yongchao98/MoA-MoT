import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a new stream to capture the output
captured_output = io.StringIO()
# Redirect stdout to the new stream
sys.stdout = captured_output

# Define the given parameters for the knowledge graph fairness audit
c = 0.85  # Confidence of the discriminatory pattern P in subgraph G
s = 0.12  # Support of the discriminatory pattern P in subgraph G
c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
beta = 0.23  # Bottleneck coefficient of the knowledge graph

# The formula to calculate the evidence ratio (E_r) is derived by comparing
# the effective evidence for the violation against the evidence suggesting it's an artifact.
#
# 1. The initial evidence for the violation is its confidence weighted by its support (c * s).
# 2. The bottleneck coefficient (beta) discounts this local evidence, as a strong bottleneck
#    (low beta) makes local patterns less representative. So, effective evidence for is (c * s * beta).
# 3. The evidence against the violation is the confidence of the counter-pattern (c_prime).
# 4. The ratio E_r is therefore (c * s * beta) / c_prime.

# Calculate the minimum ratio of evidence (E_r)
e_r = (c * s * beta) / c_prime

# --- Output Section ---
# This part of the code will print the logic, the equation, and the final result.

print("This script calculates the minimum ratio of evidence (E_r) to confirm a fairness violation.")
print("The formula considers the initial pattern's confidence and support, contradictory evidence, and the graph's structural bottleneck.")
print("\n--- Given Values ---")
print(f"Initial Confidence (c): {c}")
print(f"Initial Support (s): {s}")
print(f"Contradictory Confidence (c'): {c_prime}")
print(f"Bottleneck Coefficient (β): {beta}")

print("\n--- Calculation ---")
print("The final equation for the evidence ratio (E_r) is:")
# Output each number in the final equation
print(f"E_r = ({c} * {s} * {beta}) / {c_prime}")

# Calculate the numerator for clarity in the output
numerator = c * s * beta
print(f"E_r = {numerator:.5f} / {c_prime}")

print("\n--- Result ---")
print(f"The minimum ratio of evidence required (E_r) is: {e_r}")

# --- End of Output Section ---

# Get the captured output as a string
output_string = captured_output.getvalue()
# Restore the original stdout
sys.stdout = original_stdout
# Print the captured output
print(output_string)

# The final answer in the required format
final_answer = e_r
# <<<0.030076923076923075>>>