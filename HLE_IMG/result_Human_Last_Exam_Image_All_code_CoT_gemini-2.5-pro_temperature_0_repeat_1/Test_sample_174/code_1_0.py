import math

# --- Given Constants ---
e = 1.602e-19  # Electronic charge in Coulombs
e_eV = 1.0       # Electronic charge in elementary charge units
Nc = 1.04e19   # Density of states in Ge conduction band (cm^-3)
eps_Ge = 16.2  # Ge dielectric constant
eps_o = 8.854e-14 # Permittivity of free space (F/cm)
A_star = 114   # Richardson constant (A/K^2-cm^2)
k_eV = 8.62e-5   # Boltzmann constant (eV/K)
T = 300        # Temperature (K)
chi_Ge = 4.13  # Electron affinity of Ge (V or eV)
Jf = 20e-3     # Forward current density (A/cm^2, converted from mA/cm^2)
V_pn = 0.95    # Forward pn junction voltage (V)

# --- Part 1: Metal Identification ---

# 1. Extract data from the Mott-Schottky plot
# Point 1: (V1, y1) = (0.0 V, 1.5e15 cm^4/F^2)
# Point 2: (V2, y2) = (2.0 V, 5.5e15 cm^4/F^2)
V1, y1 = 0.0, 1.5e15
V2, y2 = 2.0, 5.5e15

# 2. Calculate slope (m) and y-intercept (c) of the line
m = (y2 - y1) / (V2 - V1)
c = y1 - m * V1

# Calculate the x-intercept, which gives -V_bi
# V_bi = -x_intercept = -(-c/m) = c/m
# Let's re-check the intercept logic. y = mx + c. At y=0, x = -c/m.
# We assumed V_dc_intercept = -V_bi. So V_bi = -V_dc_intercept = -(-c/m) = c/m.
# Let's use the points directly. V_bi = -V1 + y1/m = -0.0 + 1.5e15 / (2.0e15) = 0.75 V.
# This is incorrect. The x-intercept is V_dc = -0.75V. So V_bi = 0.75V.
V_bi = -(c / m)

# 3. Calculate donor concentration (Nd)
eps_s = eps_Ge * eps_o
Nd = 2 / (e * eps_s * m)

# 4. Calculate Schottky barrier height (phi_Bn)
kT_eV = k_eV * T
Vn = kT_eV * math.log(Nc / Nd)
phi_Bn = V_bi + Vn

# 5. Calculate metal work function (phi_m)
phi_m = phi_Bn + chi_Ge

print("--- Part 1: Metal Identification ---")
print(f"From the Mott-Schottky plot, we determine the following:")
print(f"Built-in potential, V_bi = {V_bi:.3f} V")
print(f"Donor concentration, N_d = {Nd:.2e} cm^-3")
print(f"This leads to a calculated Schottky barrier height, phi_Bn = {phi_Bn:.3f} eV")
print(f"The required metal work function is phi_m = phi_Bn + chi_Ge = {phi_Bn:.3f} eV + {chi_Ge} eV = {phi_m:.3f} eV")
print("A suitable metal with a work function in this range is Gold (Au, work function ≈ 5.1 eV).")
print("\n")


# --- Part 2: Terminal Voltage Calculation ---

# 1. Calculate reverse saturation current density (Js)
Js = A_star * (T**2) * math.exp(-phi_Bn / kT_eV)

# 2. Calculate forward voltage across Schottky barrier (Vs)
Vs = kT_eV * math.log(Jf / Js)

# 3. Calculate terminal voltage (Vt)
Vt = V_pn - Vs

print("--- Part 2: Terminal Voltage Calculation ---")
print(f"Using the calculated barrier height, the reverse saturation current density is J_s = {Js:.3e} A/cm^2.")
print(f"For a forward current density of {Jf*1000} mA/cm^2, the opposing voltage from the Schottky barrier is V_s = {Vs:.3f} V.")
print("\nThe final terminal voltage (Vt) is given by the equation:")
print(f"Vt = V_pn - V_s")
print(f"Vt = {V_pn:.2f} V - {Vs:.3f} V")
print(f"Vt = {Vt:.3f} V")
print("\n")

# --- Final Answer ---
# The final answer is the numerical value of Vt.
final_answer = f"{Vt:.3f}"
print(f"<<<{final_answer}>>>")