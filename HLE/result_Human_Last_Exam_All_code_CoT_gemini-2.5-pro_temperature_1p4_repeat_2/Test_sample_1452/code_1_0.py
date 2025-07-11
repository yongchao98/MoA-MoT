import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The name of the scientist
name = "Max Born"

# The Born rule, a fundamental postulate of quantum mechanics.
# It gives the probability density of finding a quantum particle at a given point.
# P is the probability density, and ψ (psi) is the wave function.
# The equation is P = |ψ|²
equation_latex = r"P = |\psi|^2"

# Print the output in the specified format
print(f"Name: {name}")
# The output format for the equation is LaTeX enclosed in \( ... \)
print(f"Equation: \\( {equation_latex} \\)")

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# The final answer is the generated text block
final_answer = output.strip()

# Print the final answer to the user
print(final_answer)