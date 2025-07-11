import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The user wants the equation for Max Born.
# Max Born is famous for the Born rule, which gives a probabilistic interpretation
# of the wave function. The rule states that the probability density of finding
# a particle at a point is the square of the absolute value of its wave function.
# A fundamental consequence and a common representation of this rule is the
# normalization condition: the integral of the probability density over all space is 1.

name = "Max Born"
# LaTeX representation of the normalization condition.
# This must be a string. Backslashes are escaped.
equation = "\\int |\\Psi|^2 d\\tau = 1"

# Printing the final answer in the requested format.
# This includes the numbers 2 and 1 as part of the equation string.
print(f"Name: {name}")
print(f"Equation: {equation}")

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Final output to the user
print(output)