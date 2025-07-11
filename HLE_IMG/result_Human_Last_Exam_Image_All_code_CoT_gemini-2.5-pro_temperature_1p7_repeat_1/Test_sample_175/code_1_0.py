import math

# Step 1: Define input parameters from the problem statement
# Foundation dimensions
B = 2.2  # Width of the foundation in meters
L = 2.2  # Length of the foundation in meters (B=L for square pad)
D = 2.0  # Depth of the foundation base from the ground surface in meters

# Soil properties
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2
gamma_clay_above_gwt = 18.0  # Unit weight of clay above GWT in kN/m^3
gamma_clay_below_gwt = 20.0  # Unit weight of clay below GWT in kN/m^3
gwt_depth = 0.6  # Depth of groundwater table from the surface in meters

# Partial factor for bearing resistance (EC7, DA1-1, undrained)
gamma_R_v = 1.4

print("Step 1: Input Parameters")
print(f"Foundation Width (B) = {B} m")
print(f"Foundation Length (L) = {L} m")
print(f"Foundation Depth (D) = {D} m")
print(f"Characteristic Undrained Shear Strength (c_uk) = {c_uk} kN/m^2")
print(f"Partial Factor for Bearing Resistance (γ_R,v) = {gamma_R_v}\n")

# Step 2: Calculate intermediate values
# Foundation Area (A)
A = B * L

# Total Overburden Pressure (q) at foundation base
q = (gamma_clay_above_gwt * gwt_depth) + (gamma_clay_below_gwt * (D - gwt_depth))

# Bearing Capacity Factors for undrained condition (phi_u = 0)
Nc = math.pi + 2
# Shape factor for a square foundation
s_c = 1 + 0.2 * (B / L)
# Depth factor
d_c = 1 + 0.4 * (D / B)
# Inclination factor for vertical load
i_c = 1.0

print("Step 2: Calculate intermediate values")
print(f"Foundation Area (A = B * L) = {B} * {L} = {A:.2f} m^2")
print(f"Overburden Pressure (q) = ({gamma_clay_above_gwt} * {gwt_depth}) + ({gamma_clay_below_gwt} * ({D} - {gwt_depth})) = {q:.2f} kN/m^2")
print(f"Bearing Capacity Factor (N_c = π + 2) = {Nc:.3f}")
print(f"Shape Factor (s_c = 1 + 0.2 * B/L) = 1 + 0.2 * ({B}/{L}) = {s_c:.3f}")
print(f"Depth Factor (d_c = 1 + 0.4 * D/B) = 1 + 0.4 * ({D}/{B}) = {d_c:.3f}")
print(f"Inclination Factor (i_c) = {i_c:.1f} (for vertical load)\n")

# Step 3: Calculate the characteristic ultimate bearing resistance (R_k)
# First, calculate characteristic ultimate bearing pressure (q_uk)
q_uk = (Nc * c_uk * s_c * d_c * i_c) + q
# Then, calculate characteristic bearing resistance (R_k)
R_k = q_uk * A

print("Step 3: Calculate Characteristic Bearing Resistance (R_k)")
print("The formula for characteristic ultimate bearing pressure (q_uk) is:")
print("q_uk = (N_c * c_uk * s_c * d_c * i_c) + q")
print(f"q_uk = ({Nc:.3f} * {c_uk} * {s_c:.3f} * {d_c:.3f} * {i_c:.1f}) + {q:.2f}")
q_net_term = Nc * c_uk * s_c * d_c * i_c
print(f"q_uk = {q_net_term:.2f} + {q:.2f} = {q_uk:.2f} kN/m^2")
print("\nThe formula for characteristic bearing resistance (R_k) is:")
print("R_k = q_uk * A")
print(f"R_k = {q_uk:.2f} * {A:.2f} = {R_k:.2f} kN\n")


# Step 4: Calculate the design bearing resistance (R_d)
R_d = R_k / gamma_R_v

print("Step 4: Calculate Design Bearing Resistance (R_d)")
print("The formula for design bearing resistance (R_d) is:")
print("R_d = R_k / γ_R,v")
print(f"R_d = {R_k:.2f} / {gamma_R_v}")
print(f"R_d = {R_d:.2f} kN")

print("\n--- Final Answer ---")
print("The design resistance calculation is as follows:")
print(f"R_d = (A * [ (π + 2) * c_uk * s_c * d_c * i_c + q ]) / γ_R,v")
print(f"R_d = ({A:.2f} * [({math.pi + 2:.3f}) * {c_uk} * {s_c:.3f} * {d_c:.3f} * {i_c:.1f} + {q:.2f}]) / {gamma_R_v}")
final_calculation = f"R_d = ({A:.2f} * [{q_uk:.2f}]) / {gamma_R_v}"
print(final_calculation)
final_calculation_2 = f"R_d = {R_k:.2f} / {gamma_R_v}"
print(final_calculation_2)
print(f"Final Design Resistance (R_d) = {R_d:.2f} kN")
print(f"\nFinal answer rounded to one decimal place is: {R_d:.1f} kN")