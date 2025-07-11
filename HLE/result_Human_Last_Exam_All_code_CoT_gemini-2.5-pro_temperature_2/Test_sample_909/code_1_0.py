import sympy

# Define the symbols used in the expressions
# V_0: Applied DC voltage
# r: radial distance
# pi: The mathematical constant pi
# sigma_1, sigma_2: Ohmic conductivities for region 1 and 2
# i_phi: Unit vector in the azimuthal direction in cylindrical coordinates

# Create symbolic representations for clarity
V_0, r, sigma_1, sigma_2 = sympy.symbols('V_0 r sigma_1 sigma_2')
pi = sympy.pi

# Expression for the electric field in Region 1 (0 < phi < pi/2)
# E_1 = (2 * V_0 * sigma_2 / (r * pi * (sigma_1 + sigma_2))) * i_phi
numerator_E1 = 2 * V_0 * sigma_2
denominator_E1 = r * pi * (sigma_1 + sigma_2)
E1_str = f"[ {sympy.pretty(numerator_E1)} / ({sympy.pretty(denominator_E1)}) ] * i_phi"


# Expression for the electric field in Region 2 (pi/2 < phi < pi)
# E_2 = (2 * V_0 * sigma_1 / (r * pi * (sigma_1 + sigma_2))) * i_phi
numerator_E2 = 2 * V_0 * sigma_1
denominator_E2 = r * pi * (sigma_1 + sigma_2)
E2_str = f"[ {sympy.pretty(numerator_E2)} / ({sympy.pretty(denominator_E2)}) ] * i_phi"


# Print the final results
print("The electric field in region 1 (0 < phi < pi/2) is:")
print(f"E_1 = {E1_str}\n")

print("The electric field in region 2 (pi/2 < phi < pi) is:")
print(f"E_2 = {E2_str}")