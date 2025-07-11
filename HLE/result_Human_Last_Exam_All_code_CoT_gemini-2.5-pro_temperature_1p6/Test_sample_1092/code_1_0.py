import math

# Step 1: Define baseline constants and problem variables.
# Standard JNF dose for Bluewool 3 is ~500,000 lux-hours.
jnf_dose_bw3 = 500000

# The object's Bluewool rating is 1.
object_bw_rating = 1

# The light is UV-rich. We'll use a damage factor of 2 (i.e., twice as damaging as standard).
uv_damage_factor = 2

# Exposure conditions.
lux_intensity = 50  # lux
hours_per_day = 8   # Assuming a standard 8-hour exhibition day for "daily" exposure.
days_per_year = 365.25 # To account for leap years.

# Step 2: Adjust JNF dose for the object's Bluewool rating.
# Fading rate doubles for each step down from BW 3.
# The object is (3 - 1) = 2 steps more sensitive.
sensitivity_multiplier = 2**(3 - object_bw_rating)
jnf_dose_object_standard = jnf_dose_bw3 / sensitivity_multiplier

# Step 3: Adjust for UV-rich light damage.
effective_jnf_dose = jnf_dose_object_standard / uv_damage_factor

# Step 4: Calculate daily light exposure.
daily_exposure = lux_intensity * hours_per_day

# Step 5: Calculate the time in years until one JNF.
time_in_days = effective_jnf_dose / daily_exposure
time_in_years = time_in_days / days_per_year

# Print the final equation with all the numbers.
print("This calculation determines the time in years for one Just Noticeable Fade (JNF).")
print("It is based on a series of standard values and assumptions from museum conservation science.\n")
print(f"Calculation breakdown:")
print(f"JNF Dose for Bluewool 1 (lux-hours) = {jnf_dose_bw3} / (2^(3 - {object_bw_rating})) = {jnf_dose_object_standard}")
print(f"Effective JNF Dose with UV factor (lux-hours) = {jnf_dose_object_standard} / {uv_damage_factor} = {effective_jnf_dose}")
print(f"Daily Exposure (lux-hours/day) = {lux_intensity} lux * {hours_per_day} hours = {daily_exposure}")
print("\nFinal Equation (Years):")
print(f"Years = (Effective JNF Dose / Daily Exposure) / Days per Year")
print(f"Years = ({effective_jnf_dose} / {daily_exposure}) / {days_per_year}")
print(f"Years = {time_in_days} days / {days_per_year} days_per_year")
print(f"Result: {time_in_years:.2f} years")