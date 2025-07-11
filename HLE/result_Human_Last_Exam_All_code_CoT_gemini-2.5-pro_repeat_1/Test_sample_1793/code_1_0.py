import math

# --- Constants and Inputs ---

# Beam parameters
I_current_uA = 20.0  # Proton current in microamperes
t_irr_hours = 4.0    # Irradiation time in hours
E_in_MeV = 15.0      # Initial proton energy in MeV
E_out_MeV = 12.0     # Exit proton energy in MeV

# Target parameters (natural abundance Gd2O3)
M_Gd = 157.25        # Molar mass of Gadolinium in g/mol
M_O = 15.999         # Molar mass of Oxygen in g/mol
M_Gd2O3 = 2 * M_Gd + 3 * M_O

# Product parameters
T_half_days = 5.32   # Half-life of Tb-155 in days

# Cross-section data in millibarns (mb)
cs_15MeV = 182.82
cs_14MeV = 172.16
cs_13MeV = 163.3
cs_12MeV = 150.48

# Physical constants
N_A = 6.02214076e23  # Avogadro's number (atoms/mol)
e_charge = 1.602176634e-19  # Elementary charge in Coulombs

# --- Step 1: Calculate average cross-section ---
# Average the cross-sections over the energy range
sigma_avg_mb = (cs_15MeV + cs_14MeV + cs_13MeV + cs_12MeV) / 4.0
# Convert mb to cm^2 (1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2)
sigma_avg_cm2 = sigma_avg_mb * 1e-27

# --- Step 2: Calculate target thickness in g/cm^2 ---
# Define the stopping range polynomial function
def stopping_range_g_cm2(X):
    # Y = -0.00001208736486811230*X^3 + 0.00194595770392697000*X^2 + 0.00794283377547150000*X - 0.00360695486492614000
    c3 = -0.00001208736486811230
    c2 = 0.00194595770392697000
    c1 = 0.00794283377547150000
    c0 = -0.00360695486492614000
    return c3 * X**3 + c2 * X**2 + c1 * X + c0

# Calculate range at initial and final energies
range_in = stopping_range_g_cm2(E_in_MeV)
range_out = stopping_range_g_cm2(E_out_MeV)

# The thickness is the difference in ranges
thickness_g_cm2 = range_in - range_out

# --- Step 3: Calculate the number of target atoms per cm^2 ---
# Number of Gd atoms per gram of Gd2O3 (2 Gd atoms per molecule)
num_Gd_atoms_per_gram = (2 * N_A) / M_Gd2O3
# Number of target atoms per unit area
N_target_atoms_cm2 = num_Gd_atoms_per_gram * thickness_g_cm2

# --- Step 4: Calculate the number of protons per second (flux) ---
I_current_A = I_current_uA * 1e-6  # Convert uA to A (C/s)
I_p_protons_per_sec = I_current_A / e_charge

# --- Step 5: Calculate Activity at End-of-Bombardment (EOB) ---
# Convert time units to seconds for consistency
t_irr_sec = t_irr_hours * 3600
T_half_sec = T_half_days * 24 * 3600
# Calculate decay constant (lambda)
lambda_per_sec = math.log(2) / T_half_sec

# Calculate the saturation factor
saturation_factor = 1 - math.exp(-lambda_per_sec * t_irr_sec)

# Calculate activity in Bq (disintegrations per second)
# The formula is: Activity = (Proton Flux) * (Target Atoms) * (Cross Section) * (Saturation Factor)
activity_Bq = I_p_protons_per_sec * N_target_atoms_cm2 * sigma_avg_cm2 * saturation_factor

# --- Step 6: Convert Bq to mCi ---
# 1 Ci = 3.7e10 Bq, so 1 mCi = 3.7e7 Bq
activity_mCi = activity_Bq / 3.7e7

# --- Final Output ---
print("This calculation determines the thick target yield of Tb-155.")
print("The activity (A) is calculated using the formula: A = I * N * σ * (1 - e^(-λ*t))\n")
print("The values used in the final equation are:")
print(f"I (Number of incident protons per second): {I_p_protons_per_sec:.4e} protons/s")
print(f"N (Number of target Gadolinium atoms per cm^2): {N_target_atoms_cm2:.4e} atoms/cm^2")
print(f"σ (Average cross-section): {sigma_avg_cm2:.4e} cm^2, which is {sigma_avg_mb:.2f} mb")
print(f"λ (Decay constant): {lambda_per_sec:.4e} s^-1")
print(f"t (Irradiation time): {t_irr_sec:.0f} s")
print(f"1 - e^(-λ*t) (Saturation factor): {saturation_factor:.4f}\n")
print(f"Resulting activity = ({I_p_protons_per_sec:.4e}) * ({N_target_atoms_cm2:.4e}) * ({sigma_avg_cm2:.4e}) * ({saturation_factor:.4f})")
print(f"Resulting activity = {activity_Bq:.4e} Bq\n")

print("Converting to millicuries (1 mCi = 3.7e7 Bq):")
print(f"Final calculated yield of Tb-155 = {activity_mCi:.2f} mCi")
<<<6.51>>>