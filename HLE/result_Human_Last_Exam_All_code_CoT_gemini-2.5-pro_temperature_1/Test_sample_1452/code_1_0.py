import sys
import io

# Set stdout to utf-8 to handle special characters like Ψ
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# Scientist's name
name = "Max Born"

# The famous equation associated with the scientist (Born rule)
# p = probability density, Ψ = wavefunction
# The '2' represents the squaring of the wavefunction's magnitude.
equation_latex = r"\( p = |\Psi|^2 \)"

# Print the final output
print(f"Name: {name}")
print(f"Equation: {equation_latex}")