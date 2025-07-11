import sys
import io

# Set stdout to utf-8 to handle special characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# Scientist's name
name = "Max Born"

# The famous equation associated with Max Born is the Born rule.
# It states that the probability density is the square of the magnitude of the wave function.
# The number in this equation is 2.
equation = "|ψ|²"

print(f"Name: {name}")
print(f"Equation: {equation}")