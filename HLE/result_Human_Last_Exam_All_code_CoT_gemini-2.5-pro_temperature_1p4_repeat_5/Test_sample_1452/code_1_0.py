import sys

# The unicode for the greek letter Psi (Ψ) is U+03A8
psi = '\u03A8'

# Max Born's famous equation is the Born Rule for probability density.
name = "Max Born"
# The equation is |Ψ|^2. We will represent this as a string.
equation = f"|{psi}|^2"

# Print the name and the equation.
# The number 2 is present in the final equation string.
print(f"Name: {name}")
print(f"Equation: |{psi}|²") # Using superscript two for better display
sys.stdout.flush()
