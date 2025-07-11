import math

# Step 1: Define given parameters
B = 2.2  # Foundation width in meters
L = 2.2  # Foundation length in meters
D_f = 2.0  # Foundation depth in meters
d_w = 0.6  # Groundwater level depth in meters
gamma_1 = 18.0  # Unit weight of clay above GWL in kN/m^3
gamma_2 = 20.0  # Unit weight of clay below GWL in kN/m^3
c_uk = 57.0  # Characteristic undrained shear strength in kN/m^2
# For EC7 Design Approach 1, Combination 1 (UK National Annex)
gamma_R_v = 1.0 # Partial factor for bearing resistance

# Step 2: Calculate effective foundation area (A')
A_prime = B * L

# Step 3: Calculate overburden pressure (q) at foundation base
# The pressure is from 0.6m of soil above the water table and (2.0-0.6)=1.4m of soil below
q = (d_w * gamma_1) + ((D_f - d_w) * gamma_2)

# Step 4: Determine bearing capacity factors
# Shape factor (s_c) for a square foundation (B/L = 1)
s_c = 1 + 0.2 * (B / L)

# Depth factor (d_c)
d_c = 1 + 0.4 * (D_f / B)

# Inclination factor (i_c) for vertical load
i_c = 1.0

# Step 5: Calculate ultimate bearing capacity (q_uk)
# q_uk = (pi + 2) * c_uk * s_c * d_c * i_c + q
q_uk = (math.pi + 2) * c_uk * s_c * d_c * i_c + q

# Step 6: Calculate characteristic bearing resistance (R_k)
R_k = A_prime * q_uk

# Step 7: Calculate design bearing resistance (R_d)
R_d = R_k / gamma_R_v

# Output the results
print("Calculation of the Design Bearing Resistance (R_d)\n")
print("1. Foundation and Soil Properties:")
print(f"   Foundation Width (B): {B} m")
print(f"   Foundation Length (L): {L} m")
print(f"   Foundation Depth (D_f): {D_f} m")
print(f"   Undrained Shear Strength (c_uk): {c_uk} kN/m^2")
print(f"   Partial Factor on Resistance (gamma_R;v): {gamma_R_v}\n")

print("2. Intermediate Calculations:")
print(f"   Foundation Area (A'): {A_prime:.2f} m^2")
print(f"   Overburden Pressure (q): {q:.2f} kN/m^2")
print(f"   Shape Factor (s_c): {s_c:.3f}")
print(f"   Depth Factor (d_c): {d_c:.3f}")
print(f"   Inclination Factor (i_c): {i_c:.1f}\n")

print("3. Final Calculation:")
print("   The design resistance R_d is calculated using the formula:")
print("   R_d = (A' * [ (π + 2) * c_uk * s_c * d_c * i_c + q ]) / γ_R;v\n")

print("   Substituting the values into the formula:")
# Printing the full equation with the numbers
equation = (f"   R_d = ({A_prime:.2f} * [ ({math.pi:.3f} + 2) * {c_uk} * {s_c:.3f} * "
            f"{d_c:.3f} * {i_c:.1f} + {q:.2f} ]) / {gamma_R_v:.1f}")
print(equation)

# Calculate the main term for clarity
main_term = (math.pi + 2) * c_uk * s_c * d_c * i_c
equation_step2 = (f"   R_d = ({A_prime:.2f} * [ {main_term:.2f} + {q:.2f} ]) / {gamma_R_v:.1f}")
print(equation_step2)

equation_step3 = (f"   R_d = ({A_prime:.2f} * {q_uk:.2f}) / {gamma_R_v:.1f}")
print(equation_step3)

equation_step4 = (f"   R_d = {R_k:.2f} / {gamma_R_v:.1f}")
print(equation_step4)

print("\n   The design resistance for undrained condition is:")
print(f"   R_d = {R_d:.2f} kN")
print(f"<<<{R_d:.2f}>>>")