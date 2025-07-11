import numpy as np

# --- Constants ---
e = 1.602e-19  # Electronic charge in Coulombs
eps_Ge = 16.2  # Ge dielectric constant
eps_o = 8.854e-14  # Permittivity of free space in F/cm
k_eV = 8.62e-5  # Boltzmann constant in eV/K
T = 300.0  # Temperature in Kelvin
Nc = 1.04e19  # Density of states in Ge conduction band in cm^-3
chi_Ge = 4.13  # Electron affinity of Ge in eV
A_star = 114.0  # Richardson constant in A/K^2-cm^2
V_pn = 0.95  # Forward pn junction voltage in V
Jf = 20e-3  # Forward current density in A/cm^2

# --- Part 1: Determine the Metal ---

print("--- Part 1: Finding the suitable metal ---")

# 1. Extract data from the Mott-Schottky plot
# Points from the graph: (V_dc, 1/C^2)
# Point 1 at V_dc = 0.0 V, 1/C^2 = 1.4e15 cm^4/F^2
# Point 2 at V_dc = 2.0 V, 1/C^2 = 5.6e15 cm^4/F^2
p1 = (0.0, 1.4e15)
p2 = (2.0, 5.6e15)

# Calculate slope (m) and y-intercept (y0)
m = (p2[1] - p1[1]) / (p2[0] - p1[0])
y0 = p1[1] - m * p1[0]
print(f"From the graph, the slope (m) is {m:.2e} cm^4/(F^2*V).")
print(f"The y-intercept is {y0:.2e} cm^4/F^2.")

# The x-intercept gives -Vbi. V_bi = -x_int = -(-y0/m) = y0/m
Vbi = y0 / m
# A note on convention: The plot of 1/C^2 vs V has a positive slope. Standard theory for reverse bias is 1/C^2 ~ (Vbi + Vr).
# This gives a positive slope and the x-intercept is -Vbi. We proceed with this interpretation.
print(f"The built-in potential V_bi = y0/m = {y0:.2e} / {m:.2e} = {Vbi:.3f} V.")

# 2. Calculate Donor Concentration (Nd)
eps_s = eps_Ge * eps_o
# slope m = 2 / (e * eps_s * Nd) => Nd = 2 / (e * eps_s * m)
Nd = 2 / (e * eps_s * m)
print(f"The donor concentration Nd = 2 / (e * ε_s * m) = 2 / ({e:.3e} C * {eps_s:.3e} F/cm * {m:.2e} cm^4/(F^2*V)) = {Nd:.3e} cm^-3.")

# 3. Calculate Schottky Barrier Height (phi_Bn)
kT_eV = k_eV * T
# Ec - Ef = kT * ln(Nc/Nd)
Ec_Ef = kT_eV * np.log(Nc / Nd)
print(f"The energy level difference (Ec - Ef) = kT*ln(Nc/Nd) = {kT_eV:.4f} eV * ln({Nc:.2e}/{Nd:.3e}) = {Ec_Ef:.3f} eV.")

# phi_Bn = Vbi + (Ec - Ef) in eV
phi_Bn = Vbi + Ec_Ef
print(f"The Schottky barrier height ϕ_Bn = V_bi + (Ec-Ef) = {Vbi:.3f} eV + {Ec_Ef:.3f} eV = {phi_Bn:.3f} eV.")

# 4. Calculate Metal Work Function (phi_M) and identify the metal
# phi_M = phi_Bn + chi_Ge
phi_M = phi_Bn + chi_Ge
print(f"The metal work function ϕ_M = ϕ_Bn + χ_Ge = {phi_Bn:.3f} eV + {chi_Ge:.2f} eV = {phi_M:.3f} eV.")

# Metal identification based on work function
# Common metals: Cobalt (Co) ~ 5.0 eV, Nickel (Ni) ~ 5.01-5.15 eV.
print("A metal work function of ~5.00 eV corresponds closely to Nickel (Ni) or Cobalt (Co). Nickel is a common choice for contacts on Germanium.")
print("Conclusion for Part 1: The most suitable metal is Nickel (Ni).\n")


# --- Part 2: Calculate the Terminal Voltage (Vt) ---
print("--- Part 2: Calculating the terminal voltage ---")

# 1. Calculate Saturation Current Density (Js)
# Js = A* * T^2 * exp(-phi_Bn / (kT/e))
# Note: phi_Bn and kT are in eV, so the exponent is dimensionless.
Js = A_star * T**2 * np.exp(-phi_Bn / kT_eV)
print(f"The saturation current density Js = A* * T^2 * exp(-ϕ_Bn/kT) = {A_star} * {T}^2 * exp(-{phi_Bn:.3f}/{kT_eV:.4f}) = {Js:.3e} A/cm^2.")

# 2. Calculate Schottky Diode Forward Voltage (Vsch)
# Vsch = (kT/e) * ln(Jf / Js)
Vth = kT_eV # Thermal voltage in V
V_sch = Vth * np.log(Jf / Js)
print(f"The Schottky diode voltage V_sch = Vth * ln(J_f/J_s) = {Vth:.4f} V * ln({Jf:.3f} A/cm^2 / {Js:.3e} A/cm^2) = {V_sch:.3f} V.")

# 3. Calculate Terminal Voltage (Vt)
Vt = V_pn - V_sch
print(f"The final terminal voltage Vt = V_pn - V_sch = {V_pn:.2f} V - {V_sch:.3f} V = {Vt:.3f} V.")
<<<Nickel>>>
<<<0.597>>>