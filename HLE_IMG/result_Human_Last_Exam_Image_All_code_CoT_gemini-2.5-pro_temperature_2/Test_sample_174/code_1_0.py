import math

# Define the physical constants and given parameters
e = 1.602e-19         # Electronic charge in Coulombs
Nc_Ge = 1.04e19       # Density of states in Ge conduction band in cm^-3
epsilon_Ge = 16.2     # Dielectric constant of Ge
epsilon_0 = 8.854e-14 # Permittivity of free space in F/cm
A_star = 114          # Richardson constant for Ge in A/K^2-cm^2
k_eV = 8.62e-5        # Boltzmann constant in eV/K
T = 300               # Temperature in Kelvin
chi_Ge = 4.13         # Electron affinity of Ge in V (or eV)
Jf = 20e-3            # Forward current density in A/cm^2 (20 mA/cm^2)
V_pn = 0.95           # Forward pn junction voltage in V

print("--- Part 1: Metal Identification ---")

# Step 1: Analyze the Mott-Schottky plot
# Extract two points from the linear plot of 1/C^2 vs V_dc.
# Point 1: (V1, y1) at V_dc = 0.0 V
V1, y1 = 0.0, 1.5e15  # V, cm^4/F^2
# Point 2: (V2, y2) at V_dc = 1.0 V
V2, y2 = 1.0, 3.5e15  # V, cm^4/F^2

# Calculate the slope (m) and y-intercept (b) of the line
slope = (y2 - y1) / (V2 - V1)
y_intercept = y1

# Calculate the built-in potential (V_bi) from the x-intercept (-V_bi = -b/m)
V_bi = y_intercept / slope
print("Analysis of the Mott-Schottky plot yields:")
print(f"Slope (m) = ({y2:.1e} - {y1:.1e}) / ({V2} - {V1}) = {slope:.2e} cm^4/(F^2*V)")
print(f"Built-in potential (V_bi) = {y_intercept:.1e} / {slope:.2e} = {V_bi:.2f} V")

# Step 2: Calculate donor concentration (Nd)
# From the slope equation: m = 2 / (e * epsilon_s * Nd)
epsilon_s = epsilon_Ge * epsilon_0
Nd = 2 / (e * epsilon_s * slope)
print(f"\nCalculated donor concentration (Nd):")
print(f"Nd = 2 / (e * epsilon_Ge * epsilon_0 * m) = 2 / ({e:.4e} * {epsilon_Ge} * {epsilon_0:.4e} * {slope:.2e}) = {Nd:.2e} cm^-3")

# Step 3: Calculate Fermi level position (Ec - Ef)
kT_eV = k_eV * T
Ec_minus_Ef = kT_eV * math.log(Nc_Ge / Nd)
print(f"\nCalculated Fermi level position below conduction band (Ec - Ef):")
print(f"Ec - Ef = kT * ln(Nc / Nd) = {kT_eV:.4f} eV * ln({Nc_Ge:.2e} / {Nd:.2e}) = {Ec_minus_Ef:.3f} eV")

# Step 4: Calculate the Schottky barrier height (phi_Bn)
# phi_Bn (eV) = V_bi (V) + (Ec - Ef) (eV)
phi_Bn = V_bi + Ec_minus_Ef
print(f"\nCalculated Schottky barrier height (phi_Bn):")
print(f"phi_Bn = V_bi + (Ec - Ef) = {V_bi:.2f} V + {Ec_minus_Ef:.3f} eV = {phi_Bn:.3f} eV")

# Step 5: Calculate the metal work function (phi_M) and identify the metal
# phi_M = phi_Bn + chi_Ge
phi_M = phi_Bn + chi_Ge
print(f"\nCalculated metal work function (phi_M):")
print(f"phi_M = phi_Bn + chi_Ge = {phi_Bn:.3f} eV + {chi_Ge} eV = {phi_M:.3f} eV")
print("\nConclusion for Part 1:")
print(f"The required metal must have a work function of approximately {phi_M:.3f} eV.")
print("Based on known values, Palladium (Pd, work function ~5.22 eV) is a suitable metal choice.")


print("\n\n--- Part 2: Terminal Voltage Calculation ---")

# Step 1: Calculate the reverse saturation current density (Js)
# Js = A* * T^2 * exp(-phi_Bn / kT_eV)
Js = A_star * (T**2) * math.exp(-phi_Bn / kT_eV)
print("Calculated reverse saturation current density (Js):")
print(f"Js = A* * T^2 * exp(-phi_Bn / kT) = {A_star} * {T}^2 * exp(-{phi_Bn:.3f} / {kT_eV:.4f}) = {Js:.2e} A/cm^2")

# Step 2: Calculate the forward voltage across the Schottky diode (V_Schottky)
# Jf = Js * exp(V_Schottky / V_th), where V_th = kT/e is numerically kT_eV
V_th = kT_eV
V_Schottky = V_th * math.log(Jf / Js)
print(f"\nCalculated Schottky diode forward voltage (V_Schottky) at {Jf*1000} mA/cm^2:")
print(f"V_Schottky = V_th * ln(Jf / Js) = {V_th:.4f} V * ln({Jf} / {Js:.2e}) = {V_Schottky:.3f} V")

# Step 3: Calculate the terminal voltage (Vt)
# Vt = V_pn - V_Schottky
Vt = V_pn - V_Schottky
print(f"\nCalculated final terminal voltage (Vt):")
print(f"Vt = V_pn - V_Schottky = {V_pn} V - {V_Schottky:.3f} V = {Vt:.3f} V")

# Final answers for both parts
final_answer = ('Palladium', round(Vt, 3))
print(f"\n<<<Answer Summary>>>")
print(f"Suitable Metal: {final_answer[0]}")
print(f"Terminal Voltage (Vt): {final_answer[1]} V")
<<<('Palladium', 0.334)>>>