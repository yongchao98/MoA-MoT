import math

# --- 1. Constants and Initial Parameters ---
# Beam parameters
I = 20e-6  # Beam current in Amperes (20 µA)
t_irr_hr = 4.0  # Irradiation time in hours
E_in = 15.0  # Proton entry energy in MeV
E_out = 12.0  # Proton exit energy in MeV

# Target parameters: Gadolinium(III) Oxide (Gd2O3)
M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
M_O = 15.999 # Molar mass of Oxygen in g/mol
M_Gd2O3 = 2 * M_Gd + 3 * M_O  # Molar mass of Gd2O3 in g/mol
f_Gd = 2  # Number of Gd atoms per formula unit of Gd2O3

# Product parameters (Tb-155)
t_half_days = 5.32  # Half-life of Tb-155 in days

# Physical constants
N_A = 6.02214076e23  # Avogadro's number in mol^-1
q_p = 1.60217663e-19  # Charge of a proton in Coulombs

# Conversion factors
mb_to_cm2 = 1e-27  # Millibarn to cm^2
mCi_to_Bq = 3.7e7    # millicuries to Becquerels

print("--- Step 1: Initial Parameters ---")
print(f"Proton Current (I): {I:.2e} A")
print(f"Irradiation Time (t_irr): {t_irr_hr} hours")
print(f"Half-life of Tb-155 (t_half): {t_half_days} days")
print(f"Molar Mass of Gd2O3: {M_Gd2O3:.3f} g/mol")
print("-" * 30 + "\n")


# --- 2. Calculate Decay Constant and Saturation Factor ---
# Convert times to seconds
t_irr_s = t_irr_hr * 3600
t_half_s = t_half_days * 24 * 3600

# Calculate decay constant (lambda)
lambda_val = math.log(2) / t_half_s

# Calculate saturation factor
saturation_factor = 1 - math.exp(-lambda_val * t_irr_s)

print("--- Step 2: Decay & Saturation ---")
print(f"Decay constant (λ): {lambda_val:.4e} s^-1")
print(f"Saturation factor (1 - e^(-λt)): {saturation_factor:.5f}")
print("-" * 30 + "\n")

# --- 3. Setup for Numerical Integration ---
# Stopping range equation coefficients: Y = a*X^3 + b*X^2 + c*X + d
# The derivative is dY/dX = 3a*X^2 + 2b*X + c
coeff_a = -0.00001208736486811230
coeff_b = 0.00194595770392697000
coeff_c = 0.00794283377547150000

def dY_dX(X):
    """Calculates the derivative of the stopping range polynomial dY/dX."""
    return 3 * coeff_a * X**2 + 2 * coeff_b * X + coeff_c

# Cross section data {Energy (MeV): Cross Section (mb)}
cross_sections_mb = {
    15.0: 182.82,
    14.0: 172.16,
    13.0: 163.3,
    12.0: 150.48
}
energies = sorted(cross_sections_mb.keys()) # [12.0, 13.0, 14.0, 15.0]

# --- 4. Perform Numerical Integration (Trapezoidal Rule) ---
# The integral to solve is: ∫[from E_out to E_in] (σ(E) * (dY/dX)) dE
# The unit of σ is cm^2. The unit of dY/dX is g/cm^2/MeV. The unit of dE is MeV.
# The resulting unit of the integral is grams (g).
integral_value = 0.0
for i in range(len(energies) - 1):
    E1, E2 = energies[i], energies[i+1]
    cs1_cm2, cs2_cm2 = cross_sections_mb[E1] * mb_to_cm2, cross_sections_mb[E2] * mb_to_cm2
    f1, f2 = cs1_cm2 * dY_dX(E1), cs2_cm2 * dY_dX(E2)
    integral_value += (f1 + f2) / 2.0 * (E2 - E1)

print("--- Step 3 & 4: Integral Calculation ---")
print(f"The value of ∫ σ(E) * (dY/dX) dE from {E_out} to {E_in} MeV is:")
print(f"Integral value: {integral_value:.4e} g")
print("-" * 30 + "\n")


# --- 5. Final Yield Calculation ---
# Number of protons per second
protons_per_second = I / q_p

# Number of target atoms (Gd) per gram of compound
target_atoms_per_gram = (N_A * f_Gd) / M_Gd2O3

# Full Activity Formula: A(Bq) = (I/q_p) * S * (N_A*f/M_compound) * Integral
activity_Bq = protons_per_second * saturation_factor * target_atoms_per_gram * integral_value

# Convert activity from Bq to mCi
activity_mCi = activity_Bq / mCi_to_Bq

print("--- Step 5: Final Yield (Activity) ---")
print("Activity(Bq) = (Protons/sec) * (Saturation Factor) * (Target Atoms/g) * (Integral Value)")
print("\nFinal equation with calculated numbers:")
print(f"A(Bq) = ({protons_per_second:.4e}) * ({saturation_factor:.5f}) * ({target_atoms_per_gram:.4e}) * ({integral_value:.4e})")
print(f"Calculated Activity: {activity_Bq:.4e} Bq")
print("")
print("Conversion to millicuries (mCi):")
print(f"A(mCi) = A(Bq) / (3.7e7 Bq/mCi)")
print(f"A(mCi) = {activity_Bq:.4e} / {mCi_to_Bq:.1e}")
print(f"Final Activity: {activity_mCi:.3f} mCi")
print("-" * 30)
