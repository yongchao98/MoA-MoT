import math

# Step 1: Define physical constants and given values.
# Masses and constants are rounded to three decimal places as requested.
m_delta = 1.233      # Mass of the Delta baryon in GeV
m_p = 0.938          # Mass of the proton in GeV
k = 8.617e-14        # Boltzmann's constant in GeV/K
T = 2.73             # Temperature of the CMB in K
T_base = 2.7         # Base temperature for perturbation in K
coeff_prop = 2.71    # Proportionality coefficient for the main energy term
C_pert = 1.0e-10     # Coefficient for the perturbative energy term in GeV/K^2

# Step 2: Calculate the mean energy of a CMB photon (E_gamma).
# This includes the main term and the perturbative correction.
delta_T = T - T_base
E_gamma = (coeff_prop * k * T) + (C_pert * (delta_T**2))

# Step 3: Calculate the threshold energy of the proton (E_p_thresh).
# The formula is derived from the conservation of four-momentum: s = (P_p + P_gamma)^2 = m_delta^2
# For a head-on collision at threshold, this simplifies to m_p^2 + 4*E_p*E_gamma = m_delta^2
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p_thresh = numerator / denominator

# Step 4: Output the equations with numerical values and the final result.
print("The calculation for the proton threshold energy is as follows:")

print("\n1. First, we calculate the mean photon energy (E_gamma) at T = 2.73 K.")
print("The formula for the photon energy is: E_gamma = (2.71 * k * T) + C * (T - 2.7)^2")
print(f"E_gamma = (2.71 * {k:.3e} * {T}) + {C_pert:.1e} * ({T} - {T_base})^2")
print(f"E_gamma = {E_gamma:.15e} GeV")

print("\n2. Next, we calculate the proton threshold energy (E_p) using the formula:")
print("E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma:.15e})")
print(f"E_p = ({numerator:.15f}) / ({denominator:.15e})")

print("\n3. The final result is:")
print(f"The average threshold energy for the proton is {E_p_thresh:.3e} GeV.")

# Final answer in the required format
print(f"<<<{E_p_thresh:.3e}>>>")