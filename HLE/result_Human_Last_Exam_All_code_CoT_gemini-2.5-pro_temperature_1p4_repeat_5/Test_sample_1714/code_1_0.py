import math

# This script simulates a calculation that would be performed on the hypothetical
# Wuxing decimal-based computer to find the dark matter percentage of galaxy Pandora.

# --- Step 1: Define physical constants and galaxy properties ---
# These values would be defined as 'frac' types in the Wuxing C code.
# Velocity in km/s
V = 200.0
# Radius in kiloparsecs (kpc)
R = 10.0
# Gravitational Constant in (km/s)^2 * kpc / M_sun.
# We use 4.3e-6, which can be represented on Wuxing as the fraction 43/10e-6.
G = 4.3e-6
# Luminous Mass in solar masses (M_sun), given by Luminosity (2e9) * M/L ratio (3)
M_lum = 6e9

# --- Step 2: Calculate the total mass of the galaxy ---
# The formula is M_tot = (V^2 * R) / G.
M_tot = (V**2 * R) / G

# --- Step 3: Calculate the dark matter percentage ---
# The final equation is Percentage = (1 - Luminous_Mass / Total_Mass) * 100.
# The Wuxing C program would perform these steps using its 'frac' arithmetic.
percentage_raw = (1 - M_lum / M_tot) * 100

# The result 'p' should be rounded to one decimal place.
p = round(percentage_raw, 1)

# --- Step 4: Calculate the hypothetical memory usage 'z' ---
# The Wuxing C program requires variables for: V, R, G, M_lum, M_tot, and percentage.
# This is a total of 6 variables.
# Each variable is of type 'frac', which occupies 6 decimal positions (D) of memory.
# Memory Usage (z) = 6 variables * 6 D/variable
z = 6 * 6

# --- Step 5: Output the results as requested ---
print("Calculating the Dark Matter Percentage for Pandora.")
print("The final calculation uses the formula: Percentage = (1 - M_lum / M_tot) * 100")
print("-" * 50)
print(f"Luminous Mass (M_lum) = {M_lum:.2e} M_sun")
print(f"Total Mass (M_tot)    = {M_tot:.2e} M_sun")
print(f"Dark Matter           = {p}%")
print("-" * 50)
print("The required answer in the format p:z is printed below.")
print(f"{p}:{z}")