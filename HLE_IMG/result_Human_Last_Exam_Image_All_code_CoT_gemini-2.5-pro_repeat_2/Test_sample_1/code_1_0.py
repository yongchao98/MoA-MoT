import math

# Step 1: Define constants and input values from the problem
k_B = 1.380649e-23  # Boltzmann constant in J/K
h = 6.62607015e-34   # Planck constant in J*s
R = 8.314           # Gas constant in J/(mol*K)
T = 298             # Temperature in K

# Step 2: Identify TDI, TDTS, and delta_G_r from the energy diagram
G_tdi = -152  # Gibbs energy of the TDI (P) in kJ/mol
G_tdts = 18   # Gibbs energy of the TDTS (TS-eli) in kJ/mol
delta_G_r = 10 # Overall Gibbs energy of reaction in kJ/mol

print("Step 1: Identify the rate-determining features based on the Energetic Span Model.")
print(f"The Turnover-Determining Intermediate (TDI) is P, with G = {G_tdi} kJ/mol.")
print(f"The Turnover-Determining Transition State (TDTS) is TS-eli, with G = {G_tdts} kJ/mol.")
print("The rate-determining step is governed by the energy span between P and TS-eli.")
print("-" * 30)

# Step 3: Calculate the energetic span (delta_E) in J/mol
# Formula: delta_E = G(TDTS) - G(TDI) + delta_G_r (since TDTS is before TDI)
delta_E_kJ = G_tdts - G_tdi + delta_G_r
delta_E_J = delta_E_kJ * 1000  # Convert from kJ/mol to J/mol

print("Step 2: Calculate the energetic span (δE).")
print(f"δE = G(TDTS) - G(TDI) + ΔG_r")
print(f"δE = {G_tdts} kJ/mol - ({G_tdi} kJ/mol) + {delta_G_r} kJ/mol = {delta_E_kJ} kJ/mol")
print(f"δE = {delta_E_J} J/mol")
print("-" * 30)

# Step 4: Calculate the rate constant (k) using the Eyring equation
# k = (k_B * T / h) * exp(-delta_E / (R * T))
pre_exponential_factor = (k_B * T) / h
exponent = -delta_E_J / (R * T)
k_sec = pre_exponential_factor * math.exp(exponent)

# Step 5: Convert the rate constant to hours^-1
k_hr = k_sec * 3600

print("Step 3: Calculate the rate constant (k) in hr⁻¹.")
print(f"Using the Eyring equation: k = (k_B*T/h) * exp(-δE/(R*T))")
print(f"k = {pre_exponential_factor:.3e} s⁻¹ * exp({exponent:.3f})")
print(f"k = {k_sec:.2e} s⁻¹")
print(f"k = {k_hr:.2e} hr⁻¹")
print("-" * 30)

# Final answer rounded to two significant figures
final_answer = f"{k_hr:.1e}"
print(f"The final calculated reaction rate constant to two significant figures is {k_hr:.2e} hours⁻¹.")
