import math

# Step 1: Define all given constants and physical parameters.
# All values are in GeV, K, or derived units.

# Mass of the Delta baryon in GeV
m_delta = 1.233
# Mass of the proton in GeV (rounded to three decimal places)
m_p = 0.938
# Boltzmann's constant in GeV/K (rounded to three decimal places)
k = 8.617e-14
# Specific temperature in Kelvin
T = 2.73
# Reference temperature in Kelvin for the perturbation
T_ref = 2.7
# Proportionality coefficient for the main CMB energy term
c_prop = 2.71
# Coefficient for the perturbative energy correction term in GeV*K^-2
c_pert = 1.0e-10

# Step 2: Calculate the average energy of a CMB photon (E_gamma).
# The formula is E_gamma = (c_prop * k * T) + (c_pert * (T - T_ref)^2)
delta_T = T - T_ref
E_gamma = (c_prop * k * T) + (c_pert * (delta_T**2))

# Step 3: Calculate the threshold proton energy (E_p).
# The threshold energy for p + gamma -> Delta in a head-on collision
# with a highly relativistic proton is E_p = (m_delta^2 - m_p^2) / (4 * E_gamma).
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p = numerator / denominator

# Step 4: Output the calculation and the final answer.
# As requested, output the numbers in the final equation.
print("The threshold proton energy (E_p) is calculated using the formula: E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print("\nSubstituting the values:")
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma})")
print(f"E_p = ({numerator} - {m_p**2}) / {denominator}")
print(f"E_p = {numerator} / {denominator}")


# Format the final answer in scientific notation with three decimal places.
formatted_Ep = f"{E_p:.3e}"
print(f"\nFinal calculated threshold energy: {formatted_Ep} GeV")
<<<2.201e+11>>>