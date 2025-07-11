import math

# --- Step 1: Define constants and energy values ---
# Physical constants
k_B = 1.380649e-23  # Boltzmann constant in J/K
h = 6.62607015e-34   # Planck constant in J*s
R = 8.314462       # Ideal gas constant in J/(mol*K)
T = 298.0            # Temperature in K

# Energy values from the diagram in kJ/mol
G_tdts_kj = 18.0   # G(TS-eli), the highest energy transition state
G_tdi_kj = -152.0  # G(P), the lowest energy intermediate

print("--- Analysis based on the Energetic Span Model ---")
print("1. Identifying the Turnover-Determining States:")
print(f"The highest energy transition state (TDTS) is TS-eli at {G_tdts_kj} kJ/mol.")
print(f"The lowest energy intermediate (TDI) is P at {G_tdi_kj} kJ/mol.")
print("The rate-determining factor is the energetic span between these two states.\n")

# --- Step 2: Calculate the energetic span (Gibbs free energy of activation) ---
delta_G_kj = G_tdts_kj - G_tdi_kj
delta_G_j = delta_G_kj * 1000  # Convert from kJ/mol to J/mol

print("2. Calculating the Energetic Span (ΔG‡):")
print(f"ΔG‡ = G(TDTS) - G(TDI)")
print(f"ΔG‡ = {G_tdts_kj} kJ/mol - ({G_tdi_kj} kJ/mol) = {delta_G_kj} kJ/mol")
print(f"ΔG‡ = {delta_G_j} J/mol\n")


# --- Step 3: Calculate the rate constant using the Eyring equation ---
# k = (k_B * T / h) * exp(-ΔG‡ / (R * T))
pre_exp_factor = (k_B * T) / h
exp_term = math.exp(-delta_G_j / (R * T))
k_s = pre_exp_factor * exp_term

# Convert rate constant from s^-1 to hours^-1
k_hr = k_s * 3600

print("3. Calculating the Reaction Rate Constant (k):")
print("Using the Eyring equation: k = (k_B * T / h) * exp(-ΔG‡ / (R * T))")
print("Plugging in the values:")
print(f"k = ({k_B:.4e} J/K * {T} K / {h:.4e} J*s) * exp(-{delta_G_j:.0f} J/mol / ({R:.4f} J/(mol*K) * {T} K))")
print(f"k = {k_s:.2e} s⁻¹\n")

print("4. Converting the Rate Constant to hours⁻¹:")
print(f"k (in hours⁻¹) = k (in s⁻¹) * 3600 s/hr")
print(f"k = {k_s:.2e} s⁻¹ * 3600 s/hr = {k_hr:.2e} hours⁻¹\n")

# --- Step 4: Final Answer ---
print("The rate-determining step is defined by the energetic span between TS-eli and P.")
# Format to two significant figures using an f-string
final_answer = f"{k_hr:.1e}"
print(f"The calculated reaction rate constant to two significant figures is: {final_answer} hours⁻¹")

print(f"\n<<< {final_answer} >>>")