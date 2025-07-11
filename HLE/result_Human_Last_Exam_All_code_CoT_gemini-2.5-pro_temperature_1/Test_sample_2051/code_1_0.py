import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# Define parameters from the problem for the specific case
g = 0
n_plus = 3
n_minus = 1

# --- Part (a) ---
# The property of "piecewise polynomiality" for a function representing a geometric volume
# implies that the function is continuous. The different polynomial expressions for each
# "piece" of the domain must match at the boundaries to ensure the volume does not
# jump discontinuously with an infinitesimal change in boundary lengths.
answer_a = "Yes"

# --- Part (b) ---
# The degree of the volume polynomial Z_{g,n} is related to the complex dimension
# of the corresponding moduli space of curves, which is d = 3g - 3 + n.
# The volume is a polynomial in the boundary lengths L_i, and its degree is 2*d.

# Calculate the total number of boundaries
n = n_plus + n_minus

# Calculate the degree of the polynomial
# The formula for the degree of the polynomial in terms of the boundary lengths L is 2 * (3g - 3 + n).
degree = 2 * (3 * g - 3 + n)
answer_b = degree

# --- Final Output ---
# Print the answer in the specified format
print(f"(a) {answer_a}; (b) {answer_b}")

# As requested, show the final equation for the degree calculation
print("\nThe equation for the degree calculation is:")
print(f"2 * (3 * {g} - 3 + {n}) = {degree}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)