import sys
import io

# Set stdout to utf-8 to handle special characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# The scientist's name
name = "Max Born"

# Max Born's famous equation is the Born rule, which provides the probabilistic
# interpretation of the wave function in quantum mechanics.
# The normalization condition is a direct application of this rule.
# It states that the integral of the probability density (the square of the
# wave function's magnitude) over all space must be equal to 1.
# The numbers in this equation are 2 and 1.
equation_latex = "\\int |\\Psi|^2 \\, d\\tau = 1"

# Print the name and the equation in the requested format
print(f"Name: {name}")
print(f"Equation: \\( {equation_latex} \\)")