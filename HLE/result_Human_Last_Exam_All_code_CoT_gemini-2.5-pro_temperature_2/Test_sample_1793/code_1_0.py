import math

# Step 1: Define constants and given data
# Beam and irradiation parameters
I_uA = 20.0  # Proton current in microamperes
E_in_MeV = 15.0  # Incident proton energy in MeV
E_out_MeV = 12.0  # Exit proton energy in MeV
t_irr_hr = 4.0  # Irradiation time in hours

# Product nuclide (Tb-155) properties
T_half_days = 5.32  # Half-life of Tb-155 in days

# Target (Gd2O3) properties
M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
M_O = 15.999  # Molar mass of Oxygen in g/mol

# Cross-section data in millibarns (mb)
cross_sections_mb = {
    15: 182.82,
    14: 172.16,
    13: 163.3,
    12: 150.48,
}

# Physical constants
N_A = 6.022e23  # Avogadro's number in mol^-1
e_charge = 1.602e-19  # Charge of a proton in Coulombs

# Step 2: Calculate target thickness
def calculate_range_g_cm2(X):
    """Calculates the proton range in g/cm^2 for a given energy X in MeV."""
    Y = (-0.00001208736486811230 * X**3 +
         0.00194595770392697000 * X**2 +
         0.00794283377547150000 * X -
         0.00360695486492614000)
    return Y

range_in = calculate_range_g_cm2(E_in_MeV)
range_out = calculate_range_g_cm2(E_out_MeV)
delta_x_g_cm2 = range_in - range_out

# Step 3: Calculate average cross-section
avg_cross_section_mb = sum(cross_sections_mb.values()) / len(cross_sections_mb)
# Convert cross-section from millibarns (mb) to cm^2 (1b = 1e-24 cm^2, 1mb = 1e-27 cm^2)
avg_cross_section_cm2 = avg_cross_section_mb * 1e-27

# Step 4: Calculate the number of target atoms per cm^2
M_Gd2O3 = 2 * M_Gd + 3 * M_O
# There are 2 Gd atoms per molecule of Gd2O3
n_atoms_per_molecule = 2
N_target_atoms_per_cm2 = (n_atoms_per_molecule * N_A / M_Gd2O3) * delta_x_g_cm2

# Step 5: Calculate the production rate and activity
# Convert units to be consistent (seconds)
I_A = I_uA * 1e-6  # Current in Amperes (C/s)
t_irr_s = t_irr_hr * 3600  # Irradiation time in seconds
T_half_s = T_half_days * 24 * 3600  # Half-life in seconds

# Calculate number of protons per second (flux)
n_p_per_second = I_A / e_charge

# Calculate decay constant (lambda) for Tb-155
lambda_decay = math.log(2) / T_half_s

# Calculate the production rate R (activity at saturation) in Bq
R_production_rate_Bq = n_p_per_second * N_target_atoms_per_cm2 * avg_cross_section_cm2

# Calculate the saturation factor S
saturation_factor = 1 - math.exp(-lambda_decay * t_irr_s)

# Calculate the final activity at the end of bombardment in Bq
Activity_Bq = R_production_rate_Bq * saturation_factor

# Step 6: Convert activity to millicuries (mCi)
Bq_per_mCi = 3.7e7
Activity_mCi = Activity_Bq / Bq_per_mCi

# Print the final equation and its components
print("--- Calculation of Tb-155 Yield ---")
print("\nThe final activity is calculated using the formula: Activity(mCi) = (R * S) / 3.7e7\n")
print("--- Components of the Equation ---")

print(f"1. Production Rate (R) at Saturation: R = n_p * N * <σ>")
print(f"   - Number of protons per second (n_p): {n_p_per_second:.4e} protons/s")
print(f"   - Number of target Gd atoms (N):     {N_target_atoms_per_cm2:.4e} atoms/cm^2")
print(f"   - Average Cross-Section (<σ>):       {avg_cross_section_cm2:.4e} cm^2")
print(f"   => R = {n_p_per_second:.4e} * {N_target_atoms_per_cm2:.4e} * {avg_cross_section_cm2:.4e}")
print(f"   => R = {R_production_rate_Bq:.4e} Bq\n")

print(f"2. Saturation Factor (S): S = 1 - exp(-λ * t)")
print(f"   - Decay Constant (λ):              {lambda_decay:.4e} s^-1")
print(f"   - Irradiation Time (t):            {t_irr_s:.0f} s")
print(f"   => S = 1 - exp(-{lambda_decay:.4e} * {t_irr_s:.0f})")
print(f"   => S = {saturation_factor:.4f}\n")

print("--- Final Yield Calculation ---")
print(f"Activity (mCi) = ({R_production_rate_Bq:.4e} Bq * {saturation_factor:.4f}) / {Bq_per_mCi:.1e} Bq/mCi")
print(f"Activity (Bq)  = {Activity_Bq:.4e} Bq")
print(f"Final Yield    = {Activity_mCi:.3f} mCi")

print(f"\n<<<The thick target yield of Tb-155 is {Activity_mCi:.3f} millicuries.>>>")
<<<6.497>>>