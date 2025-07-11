import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# Name of the scientist
name = "Max Born"

# The famous equation associated with Max Born (The Born Rule).
# It represents the probability density of finding a particle.
equation = "| \Psi |^2"

# Print the result in the requested format
print(f"Name: {name}")
print(f"Equation: {equation}")

# Restore stdout
sys.stdout = old_stdout
output = captured_output.getvalue()

print(output)