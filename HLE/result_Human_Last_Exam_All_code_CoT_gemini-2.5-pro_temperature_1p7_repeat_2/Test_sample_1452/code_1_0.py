import sys
import io

# Set stdout to utf-8 to handle special characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# The name of the scientist
name = "Max Born"

# The canonical commutation relation from matrix mechanics, a foundational
# part of quantum mechanics co-developed by Born.
# The format is inspired by LaTeX as seen in the examples.
equation = "[q, p] = i\\hbar"

# Print the name and the equation
print(f"Name: {name}")
print(f"Equation: \\( {equation} \\)")