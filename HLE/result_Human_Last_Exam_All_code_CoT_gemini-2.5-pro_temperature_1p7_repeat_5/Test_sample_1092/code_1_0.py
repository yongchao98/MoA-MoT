import math

# Step 1: Define the constants and variables based on conservation standards and the problem statement.
# ISO Bluewool Rating 1 shows a Just Noticeable Fade (JNF) at 100,000 lux-hours.
lux_hours_for_jnf_iso1 = 100000

# Light level the object is exposed to in lux.
lux_level = 50

# We assume a standard museum display day is 8 hours.
hours_per_day_exposure = 8

# Using 365.25 to account for leap years on average.
days_per_year = 365.25

# Step 2: Calculate the total annual light exposure in lux-hours.
annual_lux_hours = lux_level * hours_per_day_exposure * days_per_year

# Step 3: Calculate the number of years it will take to reach the JNF threshold.
years_to_fade = lux_hours_for_jnf_iso1 / annual_lux_hours

# Output the results clearly, including the equation used.
print(f"An object with ISO Bluewool Rating 1 will experience a 'Just Noticeable Fade' at {lux_hours_for_jnf_iso1:,} lux-hours.")
print(f"With a daily exposure of {lux_level} lux for {hours_per_day_exposure} hours per day:")
print(f"The number of years for the next just noticeable fade is calculated as:")
print(f"Years = {lux_hours_for_jnf_iso1} / ({lux_level} * {hours_per_day_exposure} * {days_per_year})")
print(f"Result: {years_to_fade:.2f} years")
