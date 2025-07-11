import math

# --- Step 1: Define constants and inputs ---

# Particle beam and irradiation parameters
I = 20e-6  # Current in Amperes (20 µA)
E_in = 15.0  # Initial energy in MeV
E_out = 12.0  # Exit energy in MeV
t_irr = 4.0 * 3600.0  # Irradiation time in seconds (4 hours)

# Product nuclide: Terbium-155
T_half_days = 5.32 # Half-life in days
T_half_s = T_half_days * 24.0 * 3600.0  # Half-life in seconds

# Target: Gadolinium(III) Oxide (Gd2O3)
M_Gd = 157.25  # Molar mass of Gadolinium in g/mol
M_O = 15.999   # Molar mass of Oxygen in g/mol
M_Gd2O3 = 2 * M_Gd + 3 * M_O  # Molar mass of Gd2O3 in g/mol

# Physical constants
N_A = 6.02214076e23  # Avogadro's number in mol^-1
e_charge = 1.60217663e-19  # Elementary charge in Coulombs

# Cross section data {Energy (MeV): Cross Section (mb)}
cs_data = {
    15.0: 182.82,
    14.0: 172.16,
    13.0: 163.3,
    12.0: 150.48
}

# --- Step 2: Define function for the derivative of the stopping range ---

def r_prime(E):
    """
    Calculates the derivative of the range equation (dR/dE) in g/cm^2/MeV.
    R(E) = -0.000012087...*E^3 + 0.0019459...*E^2 + 0.0079428...*E - 0.0036069...
    R'(E) = 3*(-0.000012087...)E^2 + 2*(0.0019459...)E + 0.0079428...
    """
    coeff_3 = -0.00001208736486811230
    coeff_2 = 0.00194595770392697000
    coeff_1 = 0.00794283377547150000
    return 3 * coeff_3 * E**2 + 2 * coeff_2 * E + coeff_1

# --- Step 3: Perform Numerical Integration (Trapezoidal Rule) ---
# The integral is ∫ [from E_out to E_in] of σ(E) * (dR/dE) dE

energies = sorted(cs_data.keys()) # energies = [12.0, 13.0, 14.0, 15.0]
integral_val = 0.0
for i in range(len(energies) - 1):
    E1 = energies[i]
    E2 = energies[i+1]

    # Convert cross sections from mb to cm^2 (1 mb = 1e-27 cm^2)
    sigma1_cm2 = cs_data[E1] * 1e-27
    sigma2_cm2 = cs_data[E2] * 1e-27

    # Calculate f(E) = σ(E) * R'(E)
    f1 = sigma1_cm2 * r_prime(E1)
    f2 = sigma2_cm2 * r_prime(E2)

    # Add the area of the trapezoid for the segment [E1, E2]
    segment_integral = (f1 + f2) / 2.0 * (E2 - E1)
    integral_val += segment_integral

# --- Step 4: Calculate other terms in the yield equation ---

# Number of protons per second (Ip)
I_p = I / e_charge

# Number of target Gd atoms per gram of Gd2O3 (nt)
# There are 2 Gd atoms per formula unit of Gd2O3
n_t = (2 * N_A) / M_Gd2O3

# Decay constant (lambda) for Tb-155
lambda_decay = math.log(2) / T_half_s

# Saturation factor (S)
saturation_factor = 1 - math.exp(-lambda_decay * t_irr)

# --- Step 5: Calculate the final activity ---

# Activity (A) in Bq (decays/second)
# A = Ip * nt * Integral * S
activity_Bq = I_p * n_t * integral_val * saturation_factor

# Convert Activity from Bq to mCi (1 mCi = 3.7e7 Bq)
activity_mCi = activity_Bq / 3.7e7

# --- Step 6: Print the final equation and the result ---

print("The yield (Activity) is calculated by the equation: A(Bq) = Ip * nt * Integral * S_factor")
print("Where:")
print(f"Ip (protons/sec) = {I:.1e} A / {e_charge:.4e} C = {I_p:.4e} s^-1")
print(f"nt (Gd atoms/g) = (2 * {N_A:.4e}) / {M_Gd2O3:.3f} g/mol = {n_t:.4e} atoms/g")
print(f"Integral [∫ σ(E) * (dR/dE) dE] = {integral_val:.4e} g")
print(f"S_factor [1 - e^(-λt)] = 1 - exp(-({math.log(2):.4f} / {T_half_s:.1f} s) * {t_irr:.0f} s) = {saturation_factor:.6f}")
print("\nFinal Yield Calculation (in Bq):")
print(f"A (Bq) = {I_p:.4e} * {n_t:.4e} * {integral_val:.4e} * {saturation_factor:.6f}")
print(f"A (Bq) = {activity_Bq:.4e}")

print("\nConverting to millicuries (1 mCi = 3.7e7 Bq):")
print(f"Yield (mCi) = {activity_Bq:.4e} Bq / 3.7e7 Bq/mCi")
print(f"Yield of Tb-155 = {activity_mCi:.2f} mCi")

print(f"\n<<<{activity_mCi:.2f}>>>")
