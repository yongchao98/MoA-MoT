import math

# --- User-Provided Parameters ---
E_in = 15.0  # Initial proton energy in MeV
E_out = 12.0 # Exit proton energy in MeV
proton_current_uA = 20.0 # Proton current in microAmps
irradiation_time_h = 4.0 # Total irradiation time in hours
half_life_days = 5.32 # Half-life of Tb-155 in days

# --- Material and Nuclear Data ---
# Molar masses in g/mol
molar_mass_Gd = 157.25
molar_mass_O = 15.999

# Provided cross sections in millibarns (mb)
cross_sections_mb = {
    15: 182.82,
    14: 172.16,
    13: 163.3,
    12: 150.48,
}

# --- Physical Constants ---
charge_of_proton_C = 1.602176634e-19 # Coulombs
avogadro_number = 6.02214076e23 # atoms/mol
Bq_per_mCi = 3.7e7 # Becquerels per millicurie

# --- Step-by-Step Calculation ---

# 1. Convert inputs to standard units
proton_current_A = proton_current_uA * 1e-6
irradiation_time_s = irradiation_time_h * 3600
half_life_s = half_life_days * 24 * 3600

# Calculate proton flux and decay constant
I_pps = proton_current_A / charge_of_proton_C
decay_constant_lambda = math.log(2) / half_life_s

# 2. Calculate the effective mass thickness (delta_R_mass) using the provided polynomial
def range_g_cm2(E_MeV):
    """Calculates stopping range (Y) in g/cm^2 based on proton energy (X) in MeV."""
    Y = (-0.00001208736486811230 * E_MeV**3 +
         0.00194595770392697000 * E_MeV**2 +
         0.00794283377547150000 * E_MeV -
         0.00360695486492614000)
    return Y

range_in = range_g_cm2(E_in)
range_out = range_g_cm2(E_out)
delta_R_mass = range_in - range_out

# 3. Calculate the number of target Gd atoms per cm^2
molar_mass_Gd2O3 = 2 * molar_mass_Gd + 3 * molar_mass_O
# Number of Gd atoms per gram of the Gd2O3 compound
n_atoms_per_gram = (2 * avogadro_number) / molar_mass_Gd2O3
# Total number of target atoms in the beam's path per cm^2
n_atoms_per_cm2 = n_atoms_per_gram * delta_R_mass

# 4. Calculate the average cross-section over the energy range
sigma_avg_mb = (cross_sections_mb[15] + cross_sections_mb[14] + cross_sections_mb[13] + cross_sections_mb[12]) / 4
# Convert from millibarns to cm^2 (1 mb = 1e-27 cm^2)
sigma_avg_cm2 = sigma_avg_mb * 1e-27

# 5. Calculate the production rate R (atoms/sec)
# R = (protons/sec) * (target atoms/cm^2) * (cm^2/reaction)
production_rate_R = I_pps * n_atoms_per_cm2 * sigma_avg_cm2

# 6. Calculate the saturation factor
saturation_factor = 1 - math.exp(-decay_constant_lambda * irradiation_time_s)

# 7. Calculate the final activity in Becquerels (Bq)
activity_Bq = production_rate_R * saturation_factor

# 8. Convert the final activity to millicuries (mCi)
activity_mCi = activity_Bq / Bq_per_mCi


# --- Final Output ---
print("--- Yield Calculation Breakdown ---")
print("The final activity is calculated using the formula: A = R * (1 - e^(-lambda * t)) / (Bq/mCi)")
print("where the production rate R = I_pps * n_atoms_per_cm2 * sigma_avg_cm2\n")

print(f"I_pps (proton flux) = {I_pps:.4e} p/s")
print(f"delta_R_mass (effective thickness) = {delta_R_mass:.5f} g/cm^2")
print(f"n_atoms_per_cm2 (areal density) = {n_atoms_per_cm2:.4e} atoms/cm^2")
print(f"sigma_avg_cm2 (avg. cross section) = {sigma_avg_cm2:.4e} cm^2")
print(f"lambda (decay constant) = {decay_constant_lambda:.4e} s^-1")
print(f"t (irradiation time) = {irradiation_time_s:.0f} s\n")

print("--- Equation with calculated values ---")
print(f"Production Rate (R) = {I_pps:.4e} * {n_atoms_per_cm2:.4e} * {sigma_avg_cm2:.4e} = {production_rate_R:.4e} atoms/s")
print(f"Saturation Factor = 1 - e^(-{decay_constant_lambda:.4e} * {irradiation_time_s:.0f}) = {saturation_factor:.5f}")
print(f"Activity (Bq) = {production_rate_R:.4e} * {saturation_factor:.5f} = {activity_Bq:.4e} Bq\n")

print("--- Final Tb-155 Yield ---")
print(f"The thick target yield of Tb-155 is {activity_mCi:.3f} mCi.")
