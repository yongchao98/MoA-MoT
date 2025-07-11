import math

# Step 1: Define physical constants and given values.
# As instructed, masses and physical constants are rounded to three decimal places.

# Mass of the Delta baryon in GeV.
m_delta = 1.233  

# Mass of the proton in GeV. (PDG value 0.938272... rounded to 0.938)
m_p = 0.938      

# Boltzmann's constant in GeV/K. (Value 8.617333...e-14 rounded to 8.617e-14)
k = 8.617e-14    

# Temperatures in Kelvin.
T = 2.73         
T0 = 2.7        

# Coefficients for the photon energy calculation.
c1_coeff = 2.71
c2_coeff = 1.0e-10

# Step 2: Calculate the mean energy of the CMB photon (E_gamma).
delta_T = T - T0
E_gamma = (c1_coeff * k * T) + (c2_coeff * (delta_T**2))

# Step 3: Calculate the components of the threshold energy formula.
# Numerator: m_delta^2 - m_p^2
numerator = m_delta**2 - m_p**2

# Denominator: 4 * E_gamma
denominator = 4 * E_gamma

# Step 4: Compute the threshold energy for the proton.
E_p_th = numerator / denominator

# Step 5: Print the breakdown of the calculation and the final answer.
# The instruction requires outputting each number in the final equation.
print("The threshold proton energy (E_p_th) is calculated as follows:")
print(f"E_p_th = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print(f"E_p_th = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma})")
print(f"E_p_th = ({m_delta**2} - {m_p**2}) / ({denominator})")
print(f"E_p_th = {numerator} / {denominator}")
print(f"\nThe average threshold energy for the proton is {E_p_th:.3e} GeV.")

print(f"\n<<<{E_p_th:.3e}>>>")