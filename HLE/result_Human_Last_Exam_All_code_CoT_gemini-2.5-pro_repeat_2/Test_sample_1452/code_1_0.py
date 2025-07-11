import sys
import io

# Set stdout to utf-8 to ensure Unicode characters are printed correctly
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# The scientist's name
name = "Max Born"

# The Born rule equation: ρ = |ψ|²
# ρ (rho) represents probability density.
# ψ (psi) represents the wave function.
# The equation states the probability density is the square of the magnitude of the wave function.
equation = "ρ = |ψ|²"

print(f"Name: {name}")
print(f"Equation: {equation}")