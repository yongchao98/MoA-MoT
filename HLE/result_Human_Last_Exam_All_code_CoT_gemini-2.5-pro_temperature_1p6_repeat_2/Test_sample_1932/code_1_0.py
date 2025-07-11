import sympy

# Define the symbols for the parameters
d, h, H, rho, t, g = sympy.symbols('d h H rho t g')
pi = sympy.pi

# Construct the expression for the change in weight based on the analysis
# Model: Change in weight is due to a "chain-like" fall, equal to the weight of the sand in-flight.
# Mass flow rate
m_dot = (pi * d**2 * h * rho) / (4 * t)
# Characteristic fall distance
L = H - h
# Velocity at impact
v_impact = sympy.sqrt(2 * g * L)
# Change in weight
delta_W = m_dot * v_impact

# Print the formula as a string
# The problem asks for the expression, not a numerical calculation.
# We will format the output to be readable.
final_equation = f"delta_W = (pi * d**2 * h * rho / (4 * t)) * sqrt(2 * g * (H - h))"

print("The estimated change in weight, Delta W, is given by the expression:")
print(final_equation)