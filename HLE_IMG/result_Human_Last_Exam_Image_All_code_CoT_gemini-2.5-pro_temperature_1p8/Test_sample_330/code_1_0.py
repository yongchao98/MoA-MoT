import sys
import io

# Capture original stdout to restore it later
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Analysis Logic and Answer ---

# Plot 1: Baseline simulation, showing standard coupled chaotic dynamics. Code: '0'
code_1 = "0"

# Plot 2: Denser attractor and more rapid oscillations than the baseline,
# indicating stronger chaos. This is caused by increasing the Rayleigh number.
# Code: 'R' (Rayleigh doubled)
code_2 = "R"

# Plot 3: Attractor shape is identical to the baseline (Plot 1), but the
# trajectory is different. This is the hallmark of changing the initial conditions
# in a chaotic system.
# Code: 'Z' (Z0 doubled/halved)
code_3 = "Z"

# Plot 4: The two systems (blue and orange) are almost perfectly synchronized,
# suggesting the governing equations have become symmetric. This occurs when the
# asymmetry parameter μ is set to 1.
# Code: 'M' (mu doubled, e.g., from 0.5 to 1.0)
code_4 = "M"

# Plot 5: The dynamics decay to a stable fixed point. This loss of chaos happens
# when the driving force (Rayleigh number) is reduced below a critical value.
# Code: 'r' (Rayleigh halved)
code_5 = "r"

# Plot 6: The oscillations are visibly slower than the baseline. A larger
# Prandtl number (P) slows down the evolution of the velocity variables,
# resulting in lower frequency oscillations.
# Code: 'P' (Prandtl doubled)
code_6 = "P"

# Construct the final answer string by concatenating the codes in order.
final_answer_string = code_1 + code_2 + code_3 + code_4 + code_5 + code_6

print(final_answer_string)

# --- End of analysis ---

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the buffer
output = captured_output.getvalue()

# Print the final result in the required format
print(output)
print(f"<<<{output.strip()}>>>")