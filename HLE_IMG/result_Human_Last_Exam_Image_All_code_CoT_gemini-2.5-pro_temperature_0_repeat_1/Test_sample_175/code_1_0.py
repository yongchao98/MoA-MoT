import math

# --- Given Data ---
# Soil properties
c_uk = 57  # Characteristic undrained shear strength (kN/m^2)
gamma_soil_above_wt = 18  # Unit weight of clay above water table (kN/m^3)
gamma_soil_below_wt = 20  # Unit weight of clay below water table (kN/m^3)

# Foundation geometry
B = 2.2  # Foundation width (m)
L = 2.2  # Foundation length (m)
D = 2.0  # Foundation depth (m)

# Groundwater level
d_w = 0.6  # Depth of groundwater level (m)

# Partial factor of safety (assuming a common value for undrained bearing)
gamma_R_v = 1.4

# --- Calculations ---

# Step 1: Calculate foundation area (A')
A_prime = B * L

# Step 2: Calculate total overburden pressure (q) at foundation base
# Pressure from soil above water table + pressure from soil below water table
q = (d_w * gamma_soil_above_wt) + ((D - d_w) * gamma_soil_below_wt)

# Step 3: Calculate characteristic gross bearing resistance (R_k)
# Bearing capacity factors and modifiers for undrained condition
Nc_term = math.pi + 2  # Bearing capacity term (pi + 2)
s_c = 1.2  # Shape factor for square foundation
i_c = 1.0  # Inclination factor for vertical load

# Characteristic gross bearing pressure (q_k)
q_k = (Nc_term * c_uk * s_c * i_c) + q

# Characteristic gross bearing resistance (R_k)
R_k = q_k * A_prime

# Step 4: Calculate design resistance (R_d)
R_d = R_k / gamma_R_v

# --- Output Results ---
print("Calculation of Design Bearing Resistance (Undrained Condition)")
print("-" * 60)
print("The formula for design resistance is: R_d = R_k / γ_R;v")
print("where R_k = A' * [ (π + 2) * c_uk * s_c * i_c + q ]")
print("\nIntermediate values:")
print(f"Foundation Area (A'): {A_prime:.2f} m^2")
print(f"Overburden Pressure (q): ({d_w} m * {gamma_soil_above_wt} kN/m^3) + (({D} - {d_w}) m * {gamma_soil_below_wt} kN/m^3) = {q:.2f} kN/m^2")
print(f"Characteristic Resistance (R_k): {A_prime:.2f} * [ ({math.pi:.3f} + 2) * {c_uk} * {s_c} * {i_c} + {q:.2f} ] = {R_k:.2f} kN")

print("\nFinal Equation with values:")
# Using f-string to format the final equation with all the numbers
equation_str = (f"R_d = [({A_prime:.2f}) * [ ({math.pi:.3f} + 2) * {c_uk} * {s_c} * {i_c} + {q:.2f} ]] / {gamma_R_v}")
print(equation_str)

print(f"\nDesign Resistance (R_d) = {R_d:.1f} kN")
print(f"<<<{R_d:.1f}>>>")