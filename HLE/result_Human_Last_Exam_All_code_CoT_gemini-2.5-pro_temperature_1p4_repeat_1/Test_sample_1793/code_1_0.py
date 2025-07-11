import numpy as np

# 1. Define constants and initial parameters
I_current_A = 20e-6  # Proton current in Amperes (20 µA)
E_in_MeV = 15.0      # Initial proton energy in MeV
E_out_MeV = 12.0     # Final proton energy in MeV
t_irradiation_hr = 4.0 # Irradiation time in hours
T_half_days = 5.32   # Half-life of Tb-155 in days

# Molar masses and abundances
M_Gd2O3_g_mol = 362.5  # Molar mass of Gadolinium(III) Oxide (g/mol)
f_abundance_Gd155 = 0.1480 # Natural abundance of Gd-155

# Physical constants
e_charge_C = 1.60217663e-19 # Elementary charge in Coulombs
N_A_mol = 6.02214076e23     # Avogadro's number in mol^-1
Bq_per_mCi = 3.7e7          # Becquerels per millicurie

# Cross-section data {Energy (MeV): Cross-section (mb)}
# 1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2
cross_sections_mb = {
    15: 182.82,
    14: 172.16,
    13: 163.3,
    12: 150.48
}

# Function for the derivative of the stopping range polynomial Y(X)
# dY/dX = d/dX [-0.000012087...X^3 + 0.0019459...X^2 + 0.0079428...X - 0.0036069...]
# dY/dX represents 1 / S(E), where S(E) is mass stopping power
def dYdE(E):
    """Calculates the derivative of the range polynomial at energy E (in MeV)."""
    return (-3 * 0.00001208736486811230 * E**2 +
            2 * 0.00194595770392697000 * E +
            0.00794283377547150000)

# 2. Calculate components of the yield equation

# Proton flux (particles per second)
proton_flux = I_current_A / e_charge_C

# Time conversion and decay constant (lambda)
t_irradiation_s = t_irradiation_hr * 3600
T_half_s = T_half_days * 24 * 3600
lambda_decay_constant = np.log(2) / T_half_s

# Saturation factor
saturation_factor = 1 - np.exp(-lambda_decay_constant * t_irradiation_s)

# Number of target atoms (Gd-155) per gram of Gd2O3
# There are 2 Gd atoms per molecule of Gd2O3
n_target_atoms_per_gram = (N_A_mol / M_Gd2O3_g_mol) * 2 * f_abundance_Gd155

# 3. Perform numerical integration for the integral term
energies = sorted(cross_sections_mb.keys()) # Must be sorted for integration [12, 13, 14, 15]
# Create a list of the function to be integrated: sigma(E) * dY/dE
integrand_values = []
for E in energies:
    sigma_cm2 = cross_sections_mb[E] * 1e-27  # Convert mb to cm^2
    integrand_values.append(sigma_cm2 * dYdE(E))

# Use numpy's trapezoidal rule to integrate from E_out (12) to E_in (15)
integral_val = np.trapz(y=integrand_values, x=energies)

# 4. Calculate final activity and convert units
# Activity in Bq = (protons/s) * (atoms/g) * (unitless) * (cm^2*g)
activity_Bq = proton_flux * n_target_atoms_per_gram * saturation_factor * integral_val

# Convert activity to millicuries (mCi)
activity_mCi = activity_Bq / Bq_per_mCi

# 5. Print the final results in the required format
print("Calculation of Tb-155 Yield\n")
print("The final activity is calculated using the formula:")
print("A(mCi) = [ (I_p) * (N_t) * (S) * (Integral) ] / (Bq/mCi)\n")
print("Where:")
print(f"  I_p (Proton Flux) = {proton_flux:.4e} protons/sec")
print(f"  N_t (Target Atoms/Gram) = {n_target_atoms_per_gram:.4e} atoms/gram")
print(f"  S (Saturation Factor) = {saturation_factor:.5f}")
print(f"  Integral [∫ σ(E)(dY/dE)dE] = {integral_val:.4e} cm^2*g")
print(f"  Bq/mCi Conversion Factor = {Bq_per_mCi:.1e} Bq/mCi\n")
print("Plugging in the numbers:")
print(f"A(Bq) = ({proton_flux:.4e}) * ({n_target_atoms_per_gram:.4e}) * ({saturation_factor:.5f}) * ({integral_val:.4e})")
print(f"A(Bq) = {activity_Bq:.4e} Bq\n")
print(f"A(mCi) = {activity_Bq:.4e} Bq / {Bq_per_mCi:.1e} Bq/mCi")
print(f"A(mCi) = {activity_mCi:.3f} mCi")
print("\n---")
print("Final Equation for Tb-155 Yield in millicuries:")
print(f"Yield = [ ({proton_flux:.4e}) * ({n_target_atoms_per_gram:.4e}) * ({saturation_factor:.5f}) * ({integral_val:.4e}) ] / {Bq_per_mCi:.1e} = {activity_mCi:.3f} mCi")