import math

# --- Constants ---
e = 1.602e-19  # Electronic charge (C)
eps_0 = 8.854e-14 # Permittivity of free space (F/cm)
eps_Ge = 16.2   # Germanium dielectric constant
Nc = 1.04e19    # Density of states in Ge conduction band (cm^-3)
chi_Ge = 4.13   # Electron affinity of Ge (V)
k_eV = 8.62e-5  # Boltzmann constant (eV/K)
T = 300         # Temperature (K)
A_star = 114    # Richardson constant (A/K^2-cm^2)
Jf = 0.020      # Forward current density (A/cm^2, since 20 mA/cm^2)
V_pn = 0.95     # Forward pn junction voltage (V)

# --- Part 1: Determine the Metal ---
print("--- Part 1: Identifying the Metal ---")
print("\nStep 1: Analyzing the Mott-Schottky Plot")

# Data points from the plot: (V_dc, 1/C^2)
V1, y1 = 0.0, 1.5e15
V2, y2 = 2.0, 5.5e15

# Calculate slope (m) and y-intercept (b) of the line y = mx + b
m = (y2 - y1) / (V2 - V1)
b = y1 - m * V1
print(f"The slope (m) of the 1/C^2 vs V_dc plot is {m:.2e} cm^4/(F^2*V).")
print(f"The y-intercept is {b:.2e} cm^4/F^2.")

# The x-intercept gives -V_bi. x_int = -b/m
V_x_int = -b / m
V_bi = -V_x_int
print(f"The x-intercept is V_dc = {V_x_int:.2f} V.")
print(f"Therefore, the built-in potential (V_bi) is {V_bi:.2f} V.")

print("\nStep 2: Calculating Donor Concentration (Nd)")
eps_s = eps_Ge * eps_0
# Nd = 2 / (e * eps_s * m)
Nd = 2 / (e * eps_s * m)
print(f"The donor concentration Nd is calculated as: Nd = 2 / (e * ε_s * m)")
print(f"Nd = 2 / ({e:.4e} C * {eps_s:.4e} F/cm * {m:.2e} cm^4/(F^2*V))")
print(f"Result: Nd = {Nd:.3e} cm^-3")

print("\nStep 3: Calculating Schottky Barrier Height (phi_Bn)")
# kT in eV
kT_eV = k_eV * T
# Ec_minus_Ef in eV
Ec_minus_Ef = kT_eV * math.log(Nc / Nd)
print(f"The energy difference (Ec - Ef) = kT * ln(Nc / Nd)")
print(f"(Ec - Ef) = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.3e}) = {Ec_minus_Ef:.4f} eV")

# phi_Bn in eV, V_bi in V. Numerically equivalent because of charge 'e'.
phi_Bn_eV = V_bi + Ec_minus_Ef
print(f"The Schottky barrier height (phi_Bn) = V_bi + (Ec - Ef)/e")
print(f"phi_Bn = {V_bi:.2f} V + {Ec_minus_Ef:.4f} V = {phi_Bn_eV:.4f} eV")

print("\nStep 4: Determining Metal Work Function (phi_m) and Identifying Metal")
# phi_m = phi_Bn + chi_Ge
phi_m_eV = phi_Bn_eV + chi_Ge
print(f"The metal work function (phi_m) = phi_Bn + chi_Ge")
print(f"phi_m = {phi_Bn_eV:.4f} eV + {chi_Ge:.2f} eV = {phi_m_eV:.4f} eV")
print(f"\nThe calculated metal work function is approximately {phi_m_eV:.2f} eV.")
print("This value is very close to the work function of Gold (Au, ~5.1 eV) or Palladium (Pd, ~5.12 eV).")
print("Therefore, Gold or Palladium would be a suitable metal for the diode design.")

# --- Part 2: Calculate the Terminal Voltage ---
print("\n\n--- Part 2: Calculating the Terminal Voltage (Vt) ---")

print("\nStep 1: Calculating Reverse Saturation Current Density (Js)")
# V_T = kT/e in Volts
V_T = kT_eV # Thermal voltage in V is numerically equal to kT in eV
# Js = A* * T^2 * exp(-phi_Bn / V_T)
Js = A_star * (T**2) * math.exp(-phi_Bn_eV / V_T)
print(f"The reverse saturation current density (Js) = A* * T^2 * exp(-phi_Bn / (kT/e))")
print(f"Js = {A_star} * {T}^2 * exp(-{phi_Bn_eV:.4f} / {V_T:.4f})")
print(f"Result: Js = {Js:.3e} A/cm^2")

print("\nStep 2: Calculating Schottky Diode Forward Voltage (Vs)")
# Vs = V_T * ln(Jf / Js + 1), approximated as V_T * ln(Jf / Js) since Jf >> Js
Vs = V_T * math.log(Jf / Js)
print(f"The Schottky forward voltage (Vs) = (kT/e) * ln(Jf / Js)")
print(f"Vs = {V_T:.4f} V * ln({Jf} / {Js:.3e})")
print(f"Result: Vs = {Vs:.3f} V")

print("\nStep 3: Calculating Terminal Voltage (Vt)")
# Vt = V_pn - Vs
Vt = V_pn - Vs
print(f"The terminal voltage (Vt) is given by Vt = V_pn - Vs")
print(f"Vt = {V_pn} V - {Vs:.3f} V")
print(f"Final Result: Vt = {Vt:.3f} V")

print(f"\n<<<_VALUE_>>>{Vt:.3f}")