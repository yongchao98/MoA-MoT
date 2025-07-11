import sys
import io

# Set stdout to utf-8 to handle special characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

name = "Max Born"
# The Born rule, stating the probability density is the square of the wave function's magnitude.
equation_text = "|ψ|²"

print(f"Name: {name}")
print(f"Equation: {equation_text}")