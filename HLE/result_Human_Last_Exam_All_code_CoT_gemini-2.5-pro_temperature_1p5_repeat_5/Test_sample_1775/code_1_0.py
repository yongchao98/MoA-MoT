# Define example values for the problem.
# C_tg: Top gate capacitance per unit area (e.g., in F/m^2)
# V_tg: Top gate voltage (e.g., in V)
# C_bg: Back gate capacitance per unit area (e.g., in F/m^2)
# V_bg: Back gate voltage (e.g., in V)

C_tg = 1.5e-2  # Corresponds to 1.5 uF/cm^2
V_tg = 1.0
C_bg = 0.5e-2  # Corresponds to 0.5 uF/cm^2
V_bg = -2.0

# The dielectric constant of the transistor material (eps_s) is not
# needed directly, as its effect is already accounted for in the
# capacitance per area values (C_tg and C_bg).

# Calculate the total displacement field (D) using the principle of superposition.
# The total displacement field is the sum of the fields from the top and back gates.
# D = D_tg + D_bg = C_tg * V_tg + C_bg * V_bg

D_total = C_tg * V_tg + C_bg * V_bg

# The unit of D will be in Coulombs per square meter (C/m^2)
# because (F/m^2) * V = (C/V)/m^2 * V = C/m^2.

print("The formula for the total displacement field (D) is:")
print("D = C_tg * V_tg + C_bg * V_bg")
print("\nUsing the example values, the calculation is:")
print(f"D = ({C_tg} F/m^2) * ({V_tg} V) + ({C_bg} F/m^2) * ({V_bg} V)")

D_tg_val = C_tg * V_tg
D_bg_val = C_bg * V_bg
print(f"D = ({D_tg_val} C/m^2) + ({D_bg_val} C/m^2)")
print(f"D_total = {D_total} C/m^2")
