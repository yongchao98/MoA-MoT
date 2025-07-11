import math

# Step 1: Define physical constants and problem parameters.
# As per instructions, required physical constants (proton mass, Boltzmann constant)
# are rounded to three decimal places.

# Mass of the Delta baryon (given)
m_Delta = 1.233  # GeV

# Mass of the proton, rounded to three decimal places.
m_p = 0.938  # GeV

# Boltzmann's constant, with its significand rounded to three decimal places.
k_B = 8.617e-14  # GeV/K

# Parameters for the CMB photon energy calculation
T = 2.73  # K, the temperature of the CMB
T_base = 2.7 # K, the base temperature for the perturbative correction
linear_coeff = 2.71 # Dimensionless coefficient for the linear energy term
pert_coeff = 1e-10 # GeV/K^2, coefficient for the perturbative correction term

# Step 2: Calculate the average energy of a CMB photon (E_gamma_avg).
# The formula includes a linear term and a perturbative correction.
# E_gamma_avg = (linear_coeff * k_B * T) + (pert_coeff * (T - T_base)^2)
delta_T = T - T_base
E_gamma_avg = (linear_coeff * k_B * T) + (pert_coeff * delta_T**2)

# Step 3: Calculate the proton's threshold energy (E_p) using the kinematic formula.
# E_p = (m_Delta^2 - m_p^2) / (4 * E_gamma_avg)
numerator = m_Delta**2 - m_p**2
denominator = 4 * E_gamma_avg
E_p = numerator / denominator

# Step 4: Output the details of the calculation as requested.
print("The threshold energy for the proton (E_p) is calculated using the formula:")
print("E_p = (m_Delta^2 - m_p^2) / (4 * E_gamma_avg)\n")

print("The numerical values used in this equation are:")
print(f"m_Delta = {m_Delta} GeV")
print(f"m_p = {m_p} GeV")
print(f"E_gamma_avg = {E_gamma_avg} GeV\n")

print("Substituting these values into the formula:")
print(f"E_p = ({m_Delta**2:.6f} - {m_p**2:.6f}) / (4 * {E_gamma_avg})")
print(f"E_p = {numerator:.6f} / {denominator}")
print(f"E_p = {E_p} GeV\n")

# Step 5: Format and print the final answer in scientific notation,
# rounded to three decimal places.
E_p_scientific = f"{E_p:.3e}"
print(f"The final result for the proton's threshold energy is {E_p_scientific} GeV.")

# Output the final answer in the specified format.
print(f"<<<{E_p_scientific}>>>")