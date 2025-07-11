import math

# --- User-defined variables and conservation science constants ---

# Light exposure conditions
lux_level = 50  # in lux

# Assumption for "daily" exposure in a display setting
daily_exposure_hours = 10  # hours per day
days_per_year = 365

# ISO Bluewool 1 requires ~1,000,000 lux-hours for one Just Noticeable Fade (JNF)
jnf_threshold_lux_hours = 1_000_000

# "UV-rich light" is more damaging. We use a conservative damage factor of 2.
# This means it is twice as damaging as standard low-UV light.
uv_damage_factor = 2


# --- Calculation ---

# 1. Calculate the nominal annual light dose (total lux-hours per year)
annual_lux_hours = lux_level * daily_exposure_hours * days_per_year

# 2. Calculate the effective annual light dose, including the UV damage factor
effective_annual_lux_hours = annual_lux_hours * uv_damage_factor

# 3. Calculate the number of years to reach the JNF threshold
years_to_jnf = jnf_threshold_lux_hours / effective_annual_lux_hours


# --- Output ---

print("This script calculates the time in years for a Just Noticeable Fade (JNF) to occur.")
print("The calculation is based on standard conservation principles.\n")
print("Assumptions:")
print(f"- A 'daily' exposure is {daily_exposure_hours} hours per day.")
print(f"- 'UV-rich light' has a damage factor of {uv_damage_factor} compared to standard light.")
print(f"- The threshold for a JNF for ISO Bluewool 1 is {jnf_threshold_lux_hours:,} lux-hours.\n")

print("--- Calculation Steps ---")
print("The final equation is: Years = Total Lux-Hours for JNF / (Lux * Hours/Day * Days/Year * UV Factor)\n")

# Print the final equation with all numbers substituted
print("Final Equation with values:")
print(f"Years = {jnf_threshold_lux_hours} / ({lux_level} * {daily_exposure_hours} * {days_per_year} * {uv_damage_factor})")
print(f"Years = {jnf_threshold_lux_hours} / {int(effective_annual_lux_hours)}")

# Print the final result formatted to two decimal places
print(f"\nThe time to the next just noticeable fade is: {years_to_jnf:.2f} years.")
