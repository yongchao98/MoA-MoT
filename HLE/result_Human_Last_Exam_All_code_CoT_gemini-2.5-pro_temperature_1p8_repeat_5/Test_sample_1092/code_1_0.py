import math

# Step 1: Define the constants based on conservation science standards.

# Total light exposure (dose) for a Just Noticeable Fade (JNF) for a
# highly sensitive material (ISO Bluewool 1).
jnf_dose = 120000  # units: lux-hours

# The measured light intensity.
light_level = 50  # units: lux

# The damage factor for UV-rich light compared to standard museum light.
uv_damage_factor = 2

# The exposure duration is "daily", interpreted as 8 hours per day, 365 days a year.
hours_per_day = 8
days_per_year = 365

# Step 2: Calculate the total annual light exposure.
# This accounts for the higher damage from UV-rich light.
annual_exposure = light_level * uv_damage_factor * hours_per_day * days_per_year

# Step 3: Calculate the time in years to reach the JNF dose.
years_to_fade = jnf_dose / annual_exposure

# Finally, print the equation with the values used and the result.
print("The calculation for years until the next just noticeable fade is:")
print(f"Years = Total JNF Dose / (Light Level * UV Damage Factor * Hours per Day * Days per Year)")
print(f"Years = {jnf_dose} / ({light_level} * {uv_damage_factor} * {hours_per_day} * {days_per_year})")
print(f"The result is approximately: {years_to_fade:.2f} years")