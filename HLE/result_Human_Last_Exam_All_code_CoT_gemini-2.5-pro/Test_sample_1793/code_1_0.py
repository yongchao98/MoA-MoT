import math

# Plan:
# 1. Define all the given constants and parameters.
# 2. Implement the proton range equation to find the target thickness in g/cm^2.
# 3. Calculate the average cross-section over the energy range.
# 4. Calculate the production rate (R) of Tb-155 atoms per second.
#    R = (protons/sec) * (target_atoms/cm^2) * (avg_cross_section in cm^2)
# 5. Calculate the saturation factor (S) based on the irradiation time and half-life.
#    S = 1 - exp(-lambda * t_irradiation)
# 6. Calculate the final activity at the End of Bombardment (EOB) in Becquerels (Bq).
#    Activity (Bq) = R * S
# 7. Convert the activity from Bq to millicuries (mCi) and print the result.
# 8. Print all the intermediate values used in the final calculation.

# --- Step 1: Define Constants and Parameters ---

# Beam parameters
I_current_uA = 20.0  # Proton current in microamperes
E_in = 15.0         # Initial proton energy in MeV
E_out = 12.0        # Exit proton energy in MeV
t_irradiation_hr = 4.0 # Irradiation time in hours

# Target properties
M_Gd = 157.25       # Molar mass of Gadolinium in g/mol
M_O = 16.00         # Molar mass of Oxygen in g/mol
M_Gd2O3 = 2 * M_Gd + 3 * M_O # Molar mass of Gadolinium(III) oxide

# Product properties
T_half_days = 5.32 # Half-life of Tb-155 in days

# Cross-section data in millibarns (mb)
sigma_15_mb = 182.82
sigma_14_mb = 172.16
sigma_13_mb = 163.3
sigma_12_mb = 150.48

# Physical constants
N_A = 6.022e23      # Avogadro's number (atoms/mol)
e_charge = 1.602e-19 # Elementary charge in Coulombs

# Conversion factors
Bq_per_mCi = 3.7e7  # Becquerels per millicurie

# --- Step 2: Calculate Target Thickness ---

# Equation for proton range Y (g/cm^2) as a function of energy X (MeV)
def proton_range(X):
    Y = -0.00001208736486811230 * (X**3) \
        + 0.00194595770392697000 * (X**2) \
        + 0.00794283377547150000 * X \
        - 0.00360695486492614000
    return Y

# Calculate the range at the initial and final energies
range_in = proton_range(E_in)
range_out = proton_range(E_out)

# The effective thickness of the target is the difference in ranges
target_thickness_g_cm2 = range_in - range_out

# --- Step 3: Calculate Average Cross-section ---

# Average the given cross-sections
sigma_avg_mb = (sigma_15_mb + sigma_14_mb + sigma_13_mb + sigma_12_mb) / 4

# Convert cross-section from millibarns (mb) to cm^2
# 1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2
sigma_avg_cm2 = sigma_avg_mb * 1e-27

# --- Step 4 & 5: Calculate Production Rate (R) and Saturation Factor (S) ---

# Convert units for calculation
I_current_A = I_current_uA * 1e-6 # Current in Amperes (C/s)
t_irradiation_s = t_irradiation_hr * 3600 # Time in seconds
T_half_s = T_half_days * 24 * 3600 # Half-life in seconds

# Calculate number of protons per second
protons_per_second = I_current_A / e_charge

# Calculate number of target Gd atoms per cm^2
# Number of Gd atoms per gram of Gd2O3 = (2 * N_A) / M_Gd2O3
# Number of Gd atoms per cm^2 = (Number per gram) * (thickness in g/cm^2)
n_target_atoms_per_cm2 = (2 * N_A / M_Gd2O3) * target_thickness_g_cm2

# Calculate the production rate R (atoms/sec)
R_production_rate = protons_per_second * n_target_atoms_per_cm2 * sigma_avg_cm2

# Calculate the decay constant (lambda)
lambda_decay = math.log(2) / T_half_s

# Calculate the saturation factor (S)
S_saturation_factor = 1 - math.exp(-lambda_decay * t_irradiation_s)

# --- Step 6 & 7: Calculate Final Activity ---

# Calculate activity in Becquerels (Bq) at EOB
A_Bq = R_production_rate * S_saturation_factor

# Convert activity to millicuries (mCi)
A_mCi = A_Bq / Bq_per_mCi

# --- Final Output ---

print("--- Calculation of Tb-155 Thick Target Yield ---")
print("\n--- Equation: Activity_mCi = (R * S) / Bq_per_mCi ---")
print("Where R = (protons/sec) * (target_atoms/cm^2) * (cross_section_cm^2)")
print("And S = 1 - exp(-lambda * t_irradiation)")

print("\n--- Intermediate Values ---")
print(f"Proton beam current (protons/sec): {protons_per_second:.4e}")
print(f"Irradiation time (s): {t_irradiation_s:.2f}")
print(f"Tb-155 half-life (s): {T_half_s:.2f}")
print(f"Tb-155 decay constant, lambda (s^-1): {lambda_decay:.4e}")
print(f"Target thickness (g/cm^2): {target_thickness_g_cm2:.6f}")
print(f"Number of target Gd atoms (atoms/cm^2): {n_target_atoms_per_cm2:.4e}")
print(f"Average cross-section (cm^2): {sigma_avg_cm2:.4e}")
print(f"Production rate, R (atoms/s): {R_production_rate:.4e}")
print(f"Saturation factor, S: {S_saturation_factor:.6f}")
print(f"Activity at EOB (Bq): {A_Bq:.4e}")
print(f"Conversion factor (Bq/mCi): {Bq_per_mCi:.1e}")


print("\n--- Final Calculation ---")
print(f"Activity (mCi) = ({R_production_rate:.4e} atoms/s * {S_saturation_factor:.6f}) / {Bq_per_mCi:.1e} Bq/mCi")

print(f"\nThe calculated thick target yield of Tb-155 is: {A_mCi:.3f} mCi")
<<<6.498>>>