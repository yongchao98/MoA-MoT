import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The determined sequence of parameter identifiers
# Plot 1: q_l (14)
# Plot 2: q_s (13)
# Plot 3: mu_n (3)
# Plot 4: q_f (15)
# Plot 5: c_h (12)
# Plot 6: mu_s (2)
# Plot 7: a_i (5)
# Plot 8: mu_h (8)
# Plot 9: f_s (6)
solution_sequence = "{14, 13, 3, 15, 12, 2, 5, 8, 6}"

print(solution_sequence)

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the final result in the required format
print(output.strip())
print(f'<<<{output.strip()}>>>')