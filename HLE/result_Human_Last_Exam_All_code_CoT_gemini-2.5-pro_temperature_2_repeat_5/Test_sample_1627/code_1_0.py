import sympy

# Define symbolic variables
# Pi: Internal Pressure
# a: Inner radius
Pi, a = sympy.symbols('Pi a')

# Given geometric condition
# b: Outer radius
b = 2 * a

# Lamé's formula for maximum tangential stress (sigma_t_max)
# sigma_t_max = Pi * (a**2 + b**2) / (b**2 - a**2)

# Substitute b = 2a into the formula
numerator = a**2 + b**2
denominator = b**2 - a**2

# Simplify the expressions
simplified_numerator = sympy.simplify(numerator)
simplified_denominator = sympy.simplify(denominator)

# Calculate the ratio sigma_t_max / Pi
ratio = simplified_numerator / simplified_denominator

# The result is the factor that multiplies Pi
# sigma_t_max = ratio * Pi
# So, sigma_t_max / Pi = ratio

# Print the final equation steps
print(f"The ratio of maximum tangential stress (σ_t_max) to internal pressure (P_i) is given by:")
print(f"σ_t_max / P_i = (a² + b²) / (b² - a²)")
print(f"Given that the outer radius is twice the inner radius (b = 2a), we substitute this into the equation:")
print(f"σ_t_max / P_i = (a² + (2a)²) / ((2a)² - a²)")
print(f"σ_t_max / P_i = (a² + 4a²) / (4a² - a²)")
print(f"σ_t_max / P_i = ({simplified_numerator}) / ({simplified_denominator})")
print(f"σ_t_max / P_i = {int(ratio.p)} / {int(ratio.q)}")

# Calculate and print the decimal value for comparison
decimal_value = float(ratio)
print(f"\nThe decimal value of this ratio is approximately: {decimal_value:.3f}")