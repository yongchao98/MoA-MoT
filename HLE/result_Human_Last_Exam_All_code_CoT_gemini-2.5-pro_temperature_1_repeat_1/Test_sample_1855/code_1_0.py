import math

# Step 1: Define physical constants and parameters
# Masses are in GeV. Constants are rounded to three decimal places as specified.
m_delta = 1.233  # GeV (mass of delta baryon)
m_p = 0.938      # GeV (mass of proton)
k = 8.617e-14    # GeV/K (Boltzmann's constant)

# Temperature parameters from the problem
T = 2.73         # K (target temperature)
T_base = 2.7       # K (base temperature for perturbation)

# Coefficients for the photon energy calculation
prop_coeff = 2.71     # Dimensionless proportionality coefficient
pert_coeff = 1.0e-10  # Perturbation coefficient in GeV/K^2

# Step 2: Calculate the mean energy of the CMB photon (E_gamma)
# Calculate the temperature difference for the perturbation
delta_T = T - T_base
# Calculate the total photon energy including the main term and the perturbative correction
E_gamma = (prop_coeff * k * T) + (pert_coeff * delta_T**2)

# Step 3: Use the threshold energy formula to calculate the proton energy (E_p)
# The formula is E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p_thresh = numerator / denominator

# Step 4: Print the breakdown of the calculation and the final answer
print("The formula for the proton threshold energy (E_p) is:")
print("E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)\n")
print("The values used in the equation are:")
print(f"m_delta (mass of Delta baryon) = {m_delta} GeV")
print(f"m_p (mass of proton) = {m_p} GeV")
print(f"E_gamma (energy of CMB photon) = {E_gamma:.6e} GeV\n")

print("Substituting these values into the equation:")
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma:.6e})")
print(f"E_p = ({m_delta**2:.6f} - {m_p**2:.6f}) / ({denominator:.6e})")
print(f"E_p = {numerator:.6f} / {denominator:.6e}\n")

# Format the final result to scientific notation with three decimal places
final_answer_formatted = f"{E_p_thresh:.3e}"
print(f"The calculated average threshold energy for the proton is {final_answer_formatted} GeV.")

# Output the final answer in the specified format
print(f"<<<{final_answer_formatted}>>>")