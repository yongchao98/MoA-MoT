# Plan:
# 1. Define the standard total light exposure (in lux-hours) required to cause a "just noticeable fade" (JNF) for a material with an ISO Bluewool rating of 1.
# 2. Define the given daily light exposure conditions: intensity (lux) and duration (hours per day).
# 3. Calculate the annual light exposure based on the daily conditions.
# 4. Divide the total JNF exposure by the annual exposure to find the number of years until fading occurs.
# 5. Print the full equation and the final result.

# --- Parameters ---

# Total cumulative exposure in lux-hours to cause a just noticeable fade for Bluewool 1 material.
total_lux_hours_for_jnf = 500000

# Given light intensity in lux.
lux_intensity = 50

# Assumed hours of exposure per day (standard for a museum/gallery).
hours_per_day = 8

# Number of days in a year.
days_per_year = 365

# --- Calculation ---

# Calculate the annual exposure.
annual_exposure = lux_intensity * hours_per_day * days_per_year

# Calculate the number of years until the next just noticeable fade.
years_to_fade = total_lux_hours_for_jnf / annual_exposure

# --- Output ---

# Print the final equation with all the values and the result.
# The result is rounded to two decimal places for clarity.
print(f"To calculate the years until a just noticeable fade occurs, we use the formula:")
print(f"Total JNF Exposure / (Lux Intensity * Hours per Day * Days per Year)")
print(f"The calculation is:")
print(f"{total_lux_hours_for_jnf} / ({lux_intensity} * {hours_per_day} * {days_per_year}) = {years_to_fade:.2f} years")
