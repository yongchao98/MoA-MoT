import math

# Define constants
k_B = 1.380649e-23  # Boltzmann constant in J/K
h = 6.62607015e-34   # Planck constant in J*s
R = 8.3144626      # Ideal gas constant in J/(mol*K)
T = 298.0            # Temperature in Kelvin

# Step 1: Identify Rate-Determining Step and Energetic Span
# Gibbs free energies from the diagram in kJ/mol
intermediates = {'Cat': 0, 'Cat-GB': -61, 'R': -55, 'P5': -22}
transition_states = {'TS-add': 10, 'TS-eli': 18}

# Identify the Turnover Determining Transition State (TDTS) as the highest-energy TS
G_TDTS = transition_states['TS-eli']
TDTS_name = 'TS-eli'

# Identify the Turnover Determining Intermediate (TDI) as the lowest-energy intermediate before the TDTS
# The intermediates appearing before TS-eli are Cat, Cat-GB, R, and P5.
precursors_to_TDTS = {'Cat': 0, 'Cat-GB': -61, 'R': -55, 'P5': -22}
G_TDI = min(precursors_to_TDTS.values())
TDI_name = [name for name, G in precursors_to_TDTS.items() if G == G_TDI][0]

# The rate-determining step is the one involving the TDTS
RDS_description = f"the elimination step involving {TDTS_name} (from P5 to P)"

# Calculate the energetic span (δE) in kJ/mol and convert to J/mol
delta_E_kJ = G_TDTS - G_TDI
delta_E_J = delta_E_kJ * 1000

print("Step 1: Identify the Rate-Determining Step (RDS) and Energetic Span (δE)")
print("-------------------------------------------------------------------------")
print(f"The Turnover Determining Transition State (TDTS) is {TDTS_name} at {G_TDTS} kJ/mol.")
print(f"The Turnover Determining Intermediate (TDI) is {TDI_name} at {G_TDI} kJ/mol.")
print(f"Based on the Energetic Span Model, the rate-determining step is {RDS_description}.")
print("\nThe energetic span (δE) is the energy difference between the TDTS and the TDI:")
print(f"δE = G({TDTS_name}) - G({TDI_name})")
print(f"δE = {G_TDTS} kJ/mol - ({G_TDI} kJ/mol) = {delta_E_kJ} kJ/mol")
print("\n")

# Step 2: Calculate the reaction rate constant (k)
print("Step 2: Calculate the reaction rate constant (k) using the Eyring equation")
print("-------------------------------------------------------------------------")
print("k = (k_B * T / h) * exp(-δE / (R * T))")

# Calculate the rate constant in s^-1
pre_exp_factor = (k_B * T) / h
exponent = -delta_E_J / (R * T)
k_s = pre_exp_factor * math.exp(exponent)

print("\nPlugging in the values:")
print(f"k = (({k_B:.4e} J/K * {T} K) / {h:.4e} J*s) * exp(-{delta_E_J} J/mol / ({R:.4f} J/(mol*K) * {T} K))")
print(f"k = ({pre_exp_factor:.4g} s⁻¹) * exp({exponent:.4f})")
print(f"k ≈ {k_s:.4g} s⁻¹")
print("\n")

# Step 3: Convert to hours^-1 and round to two significant figures
k_hr = k_s * 3600
k_hr_rounded = float(f"{k_hr:.2g}")

print("Step 3: Convert to final units and round")
print("-----------------------------------------")
print(f"Converting from s⁻¹ to hr⁻¹: k (hr⁻¹) = k (s⁻¹) * 3600 s/hr")
print(f"k ≈ {k_s:.4g} s⁻¹ * 3600 s/hr ≈ {k_hr:.4g} hr⁻¹")
print(f"\nThe final reaction rate constant rounded to two significant figures is {k_hr_rounded} hr⁻¹.")
