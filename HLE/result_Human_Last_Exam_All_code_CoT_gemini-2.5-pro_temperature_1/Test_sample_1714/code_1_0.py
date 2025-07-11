import math

# Pandora galaxy parameters
velocity_kms = 200  # Velocity in km/s
radius_kpc = 10  # Radius in kpc
luminosity_Lsun = 2e9  # Luminosity in solar units
mass_to_light_ratio = 3  # Mass-to-light ratio in solar units

# Gravitational constant in units of (kpc * (km/s)^2 / M_sun)
G = 4.30091e-6

# --- Calculations ---

# 1. Calculate the total mass of the galaxy within the given radius
# Formula: M_total = (velocity^2 * radius) / G
total_mass = (velocity_kms**2 * radius_kpc) / G

# 2. Calculate the luminous mass (mass of stars, etc.)
# Formula: M_luminous = Luminosity * Mass-to-Light ratio
luminous_mass = luminosity_Lsun * mass_to_light_ratio

# 3. Calculate the dark matter mass
dark_mass = total_mass - luminous_mass

# 4. Calculate the percentage of dark matter
dark_matter_percentage = (dark_mass / total_mass) * 100

# --- Output ---
# Print the final equation with values, as requested.
# The numbers are formatted in scientific notation for clarity.
print("Equation: Dark Matter % = ((Total Mass - Luminous Mass) / Total Mass) * 100")
print(f"Calculation: (({total_mass:.2e} - {luminous_mass:.2e}) / {total_mass:.2e}) * 100 = {dark_matter_percentage:.1f}%")

# The final answer p:z requires the percentage (p) and memory usage (z).
# p is the calculated percentage, rounded to 0.1 decimal.
p = round(dark_matter_percentage, 1)

# z is the memory usage in 'D' for the hypothetical C program on Wuxing.
# The C program would require variables for the calculation:
# 5 input fracs (v, r, G, L, M_L_ratio): 5 * 6D = 30D
# 4 intermediate/result fracs (M_total, M_lum, ratio, percentage): 4 * 6D = 24D
# For printing rounded to 0.1%, we need one integer variable to hold the result * 10.
# 1 int (for printing): 1 * 5D = 5D
# Total memory z = 30D + 24D + 5D = 59D.
z = 59

# The final answer will be provided in the required format below.
# print(f"Final Answer Format: {p}:{z}")