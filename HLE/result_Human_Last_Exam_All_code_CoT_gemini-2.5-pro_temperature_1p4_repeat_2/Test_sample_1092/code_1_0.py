import math

# --- Constants and Assumptions based on Conservation Science ---

# The standard light exposure dose (in lux-hours) to cause a
# Just Noticeable Fade (JNF) on a Bluewool Scale #1 sample.
JNF_THRESHOLD_BLUEWOOL_1 = 50000

# An assumed damage multiplier for "UV-rich" light compared to standard
# UV-filtered museum lighting. Fading will occur this many times faster.
UV_DAMAGE_FACTOR = 3

# An assumed number of hours per day for "daily" exposure,
# typical for a museum or gallery setting.
HOURS_PER_DAY = 8

# Days in a year for the final calculation.
DAYS_PER_YEAR = 365

# --- Given information from the user ---
lux_level = 50

# --- Calculation ---

# 1. Calculate the effective tolerance of the object in the specified UV-rich light.
effective_tolerance_lux_hours = JNF_THRESHOLD_BLUEWOOL_1 / UV_DAMAGE_FACTOR

# 2. Calculate the total annual light exposure (dose) the object receives.
annual_exposure_lux_hours = lux_level * HOURS_PER_DAY * DAYS_PER_YEAR

# 3. Calculate the time in years until a JNF occurs.
years_to_jnf = effective_tolerance_lux_hours / annual_exposure_lux_hours

# --- Output the result ---

print("To find the time to the next Just Noticeable Fade (JNF), we follow these steps:")
print("1. Determine the JNF threshold for the material (Bluewool 1).")
print("2. Adjust the threshold for the damaging effect of UV-rich light.")
print("3. Calculate the total annual light exposure based on the given conditions.")
print("4. Divide the adjusted threshold by the annual exposure.")
print("\n--- Calculation Details ---")
print(f"Standard JNF Threshold for Bluewool 1: {JNF_THRESHOLD_BLUEWOOL_1} lux-hours")
print(f"Assumed UV Damage Factor: {UV_DAMAGE_FACTOR}")
print(f"Assumed Daily Exposure: {HOURS_PER_DAY} hours")
print(f"Light Level: {lux_level} lux")
print("\n--- Final Equation ---")
# The final requested format is to show the equation with numbers
print(f"Years = (JNF Threshold / UV Factor) / (Lux Level * Hours per Day * Days per Year)")
print(f"Years = ({JNF_THRESHOLD_BLUEWOOL_1} / {UV_DAMAGE_FACTOR}) / ({lux_level} * {HOURS_PER_DAY} * {DAYS_PER_YEAR})")
print(f"Years = {effective_tolerance_lux_hours:.2f} / {annual_exposure_lux_hours}")
print(f"Result: {years_to_jnf:.3f} years")
print("\nBased on these calculations, the next just noticeable fade will occur in approximately 0.114 years.")
