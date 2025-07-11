import math

# Plan:
# 1. Define the physical constants and given values.
# 2. Use a version of the gravitational constant G with units appropriate for the problem (kpc, km/s, solar mass).
#    This simplifies the calculation and makes it feasible for the Wuxing 'frac' type.
# 3. Calculate total mass using M_total = (v^2 * r) / G.
# 4. Calculate luminous mass using M_lum = (Mass/Light Ratio) * Luminosity.
# 5. Calculate the dark matter percentage.
# 6. Calculate the memory footprint of the equivalent Wuxing C program.
# 7. Print the results as requested.

# --- Step 1 & 2: Define constants and givens ---
# Velocity curve at radius r
v_kms = 200  # km/s

# Radius
r_kpc = 10  # kpc

# Gravitational Constant in convenient units for this problem.
# G = 4.302e-6 kpc*(km/s)^2 / M_sun
# In the Wuxing C program, this would be represented as a `frac`,
# e.g., G = 43/1e-7 {n=43, d=1, e=-7}
G = 4.302e-6

# Given luminosity and mass-to-light ratio
luminosity_solar = 2e9  # L_sun
mass_to_light_solar = 3  # M_sun/L_sun

# --- Step 3 & 4: Mass Calculations ---
# Luminous Mass in Solar Masses (M_sun)
m_lum_solar = luminosity_solar * mass_to_light_solar

# Total Mass in Solar Masses (M_sun)
m_total_solar = (v_kms**2 * r_kpc) / G

# --- Step 5: Dark Matter Percentage ---
# Dark Mass in Solar Masses
m_dark_solar = m_total_solar - m_lum_solar

# Dark Matter Percentage
dark_matter_percentage = (m_dark_solar / m_total_solar) * 100

# --- Step 6: Wuxing Memory Calculation ---
# Variables needed for the C program:
# 4 `frac` types: M_total, M_lum, G, one for calculation/result.
# sizeof(frac) = 2D (char) + 2D (char) + 2D (char) = 6D
mem_frac_vars = 4 * 6 # in D

# 3 `int` types for printing logic (val_times_10, int_part, frac_part)
# sizeof(int) = 5D
mem_int_vars = 3 * 5 # in D

total_memory_usage_D = mem_frac_vars + mem_int_vars

# --- Step 7: Print the final output ---
# Per instructions, the final code needs to output each number in the final equation.
print("--- Calculation Breakdown ---")
print("Equation: Dark Matter % = ((M_total - M_lum) / M_total) * 100")
print(f"where M_total = (v^2 * r) / G")
print(f"v (velocity) = {v_kms} km/s")
print(f"r (radius) = {r_kpc} kpc")
print(f"G (gravitational constant) = {G:.4e} kpc*(km/s)^2/M_sun")
print(f"M_lum (luminous mass) = {m_lum_solar:.1e} Solar Masses")
print(f"Resulting M_total (total mass) = {m_total_solar:.3e} Solar Masses")
print(f"Resulting M_dark (dark mass) = {m_dark_solar:.3e} Solar Masses")
print("")
print("--- Final Answer ---")

# Format the final answer as p:z
p_formatted = f"{dark_matter_percentage:.1f}%"
z_formatted = f"{total_memory_usage_D}D"

print(f"The percentage of dark matter is {p_formatted} and the required memory is {z_formatted}.")
# The final tagged answer should strictly follow the p:z format.
print(f"Final formatted answer: {p_formatted}:{z_formatted}")