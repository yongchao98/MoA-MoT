import math

# Step 1: Define the initial parameters from the problem description.
# We will use integers to represent the numbers to align with the Wuxing's decimal, non-floating-point architecture.

# Luminosity is 2e9 times the Sun's luminosity
luminosity_factor = 2 * (10**9)

# Mass/Light ratio is 3 times the Sun's ratio
mass_to_light_ratio = 3

# Radius is 10 kpc
r_kpc = 10

# Velocity is 200 km/s
v_kms = 200

# The constant for total mass calculation (2.32e5) is represented as an integer 232000
mass_formula_const = 232000

# Step 2: Calculate the Luminous Mass in solar masses.
# M_lum = luminosity * mass_to_light_ratio
m_lum_solar = luminosity_factor * mass_to_light_ratio

# Step 3: Calculate the Total Mass in solar masses using the provided formula.
# M_total = const * r * v^2
m_total_solar = mass_formula_const * r_kpc * (v_kms ** 2)

# Step 4: Calculate the Dark Mass.
m_dark_solar = m_total_solar - m_lum_solar

# Step 5: Calculate the percentage of dark matter.
percentage_dark_matter = (m_dark_solar / m_total_solar) * 100

# Round the percentage to one decimal place as required.
p = round(percentage_dark_matter, 1)

# Step 6: Calculate the memory usage (z).
# A 'frac' type uses 3 * 'char' (2D each) = 6D.
# We need variables for M_lum, M_total, M_dark, and the final percentage.
# Thus, we need 4 'frac' variables.
num_variables = 4
mem_per_variable_D = 6
z = num_variables * mem_per_variable_D

# Step 7: Print the final output.
# The prompt requires printing the numbers in the final equation.
print("Calculation Steps:")
print(f"Luminous Mass = {luminosity_factor:.0e} L_sun * {mass_to_light_ratio} M_sun/L_sun = {m_lum_solar:.2e} M_sun")
print(f"Total Mass = {mass_formula_const} * {r_kpc} kpc * ({v_kms} km/s)^2 = {m_total_solar:.2e} M_sun")
print("\nDark Matter Percentage Equation:")
print(f"( {m_total_solar:.2e} - {m_lum_solar:.2e} ) / {m_total_solar:.2e} * 100 = {p}%")

print("\nFinal Answer (p:z):")
print(f"{p}:{z}")
