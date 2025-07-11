import math

# Step 1: Define Constants and Inputs
# -------------------------------------
# Beam properties
proton_current_A = 20e-6  # 20 µA in Amperes (C/s)
E_in_MeV = 15.0         # Initial proton energy in MeV
E_out_MeV = 12.0        # Exit proton energy in MeV

# Target properties
target_density_g_cm3 = 7.41 # Density of Gd2O3 in g/cm^3
# Molar mass of Gd2O3 (Gd=157.25 g/mol, O=16.00 g/mol)
M_Gd2O3_g_mol = 2 * 157.25 + 3 * 16.00
# Number of Gd atoms per molecule of Gd2O3
k_atoms_per_molecule = 2

# Irradiation parameters
irradiation_time_hr = 4.0
T_half_life_days = 5.32

# Cross-section data in millibarns (mb)
cross_sections_mb = {
    12: 150.48,
    13: 163.3,
    14: 172.16,
    15: 182.82,
}

# Physical constants
e_charge_C = 1.60217663e-19  # Elementary charge in Coulombs
N_A_mol = 6.02214076e23     # Avogadro's number in atoms/mol
Bq_per_mCi = 3.7e7          # Conversion factor for Bq to mCi


# Step 2: Calculate Derived Constants
# -----------------------------------
# Convert times to seconds for consistency
irradiation_time_s = irradiation_time_hr * 3600
T_half_life_s = T_half_life_days * 24 * 3600

# Proton flux (particles per second)
proton_flux_per_s = proton_current_A / e_charge_C

# Decay constant (lambda) in s^-1
decay_constant_s = math.log(2) / T_half_life_s

# Saturation factor (dimensionless)
saturation_factor = 1 - math.exp(-decay_constant_s * irradiation_time_s)

# Number of target (Gd) atoms per gram of compound
n_atoms_per_gram = (k_atoms_per_molecule * N_A_mol) / M_Gd2O3_g_mol


# Step 3: Implement Stopping Power Calculation
# --------------------------------------------
# The derivative of the range polynomial Y(X) gives dR/dE
# Y = -0.000012087*X^3 + 0.0019460*X^2 + 0.0079428*X - 0.0036070
# dY/dX = dR/dE = -3*0.000012087*X^2 + 2*0.0019460*X + 0.0079428
def calculate_dR_dE(energy_MeV):
    """Calculates dR/dE in (g/cm^2)/MeV from the given polynomial derivative."""
    c3 = -0.00001208736486811230
    c2 = 0.00194595770392697000
    c1 = 0.00794283377547150000
    return 3 * c3 * energy_MeV**2 + 2 * c2 * energy_MeV + c1


# Step 4: Perform Numerical Integration (Trapezoidal Rule)
# --------------------------------------------------------
# We want to calculate Integral[ sigma(E) / S(E) dE ]
# which is equivalent to Integral[ sigma(E) * (dR/dE) dE ]
integral_value = 0.0
energy_points = sorted(cross_sections_mb.keys()) # [12, 13, 14, 15]

# The integral is approximated by summing up the areas of trapezoids
# for each 1 MeV energy step from 12 to 15 MeV.
for i in range(len(energy_points) - 1):
    E1 = energy_points[i]
    E2 = energy_points[i+1]
    
    # Get cross-sections for the two energy points and convert from mb to cm^2
    sigma1_cm2 = cross_sections_mb[E1] * 1e-27
    sigma2_cm2 = cross_sections_mb[E2] * 1e-27

    # Calculate the integrand f(E) = sigma(E) * dR/dE(E) at each point
    f1 = sigma1_cm2 * calculate_dR_dE(E1)
    f2 = sigma2_cm2 * calculate_dR_dE(E2)

    # Area of one trapezoid = (f1 + f2) / 2 * (E2 - E1)
    h = E2 - E1
    integral_value += (f1 + f2) / 2.0 * h


# Step 5: Calculate Final Activity
# --------------------------------
# Production rate at saturation (R_prod) = Flux * n_atoms_per_gram * Integral
# Activity (A) = R_prod * Saturation_Factor
activity_Bq = proton_flux_per_s * n_atoms_per_gram * integral_value * saturation_factor


# Step 6: Convert to Millicuries and Display
# ------------------------------------------
activity_mCi = activity_Bq / Bq_per_mCi

print("Calculation of Tb-155 Thick Target Yield\n")
print("The final activity is calculated using the formula:")
print("A(mCi) = (Flux * N_atoms_per_gram * Integral * Saturation_Factor) / Bq_per_mCi\n")
print("Where the values are:")
print(f"Proton Flux             = {proton_flux_per_s:.4e} protons/s")
print(f"Target Atoms per Gram   = {n_atoms_per_gram:.4e} atoms/g")
print(f"Saturation Factor       = {saturation_factor:.4f}")
print(f"Integral(sigma*dR/dE)   = {integral_value:.4e} g")
print(f"Conversion Factor       = {Bq_per_mCi:.4e} Bq/mCi\n")
print("Plugging in the numbers:\n")

# Display the final equation with all numbers plugged in
equation_str = (
    f"A(mCi) = (({proton_flux_per_s:.4e}) * ({n_atoms_per_gram:.4e}) * "
    f"({integral_value:.4e}) * ({saturation_factor:.4f})) / ({Bq_per_mCi:.4e})"
)
print(equation_str)

final_yield_str = f"Final calculated yield = {activity_mCi:.2f} mCi"
print("\n" + "="*len(final_yield_str))
print(final_yield_str)
print("="*len(final_yield_str))

print(f"\n<<<{activity_mCi:.2f}>>>")