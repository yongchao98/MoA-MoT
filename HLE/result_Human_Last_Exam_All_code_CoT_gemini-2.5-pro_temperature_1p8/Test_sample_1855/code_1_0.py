import math

# This script calculates the threshold energy for a proton-photon reaction.

# Step 1: Define physical constants with the specified precision.
m_delta = 1.233  # Mass of the Delta baryon in GeV
m_p = 0.938      # Mass of the proton in GeV (rounded to 3 decimal places from 938.272 MeV)
T = 2.73         # Given CMB temperature in Kelvin
T0 = 2.7         # Base temperature for the perturbative correction in Kelvin
k = 8.617e-14    # Boltzmann constant in GeV/K (rounded to 3 decimal places from 8.617333e-14)
C1_coeff = 2.71  # Proportionality coefficient for the main energy term
C2 = 1.0e-10     # Perturbative coefficient in GeV/K^2

# Step 2: Calculate the mean energy of a CMB photon at T=2.73 K.
# The formula is E_gamma = (2.71 * k * T) + C2 * (T - T0)^2.
delta_T = T - T0
main_energy_term = C1_coeff * k * T
perturbative_term = C2 * delta_T**2
E_gamma = main_energy_term + perturbative_term

# Step 3: Calculate the proton threshold energy (E_th).
# For a head-on collision, the threshold energy E_th can be found using the formula
# derived from the invariance of the four-momentum squared, which gives:
# E_th ≈ (m_delta^2 - m_p^2) / (4 * E_gamma).
# This approximation is extremely accurate for this case.
m_delta_sq = m_delta**2
m_p_sq = m_p**2
E_th = (m_delta_sq - m_p_sq) / (4 * E_gamma)

# Step 4: Print the equations with all numerical values to show the calculation steps.
print("First, the average CMB photon energy (E_gamma) at T=2.73 K is calculated:")
print(f"E_gamma = ({C1_coeff} * {k:.3e} * {T}) + ({C2:.1e} * ({T} - {T0})^2) = {E_gamma:.6e} GeV")

print("\nThen, this energy is used to find the proton threshold energy (E_th):")
# This line fulfills the requirement to output each number in the final equation.
print(f"E_th = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma:.6e}) = {E_th:.3e} GeV")
