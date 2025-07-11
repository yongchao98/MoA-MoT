import math

# --- Step 1: Define the constants based on conservation standards ---
# Total lux-hours of exposure to cause a "just noticeable fade" (JNF) for an
# object with an ISO Bluewool rating of 1.
total_lux_hours_for_jnf = 50000

# Exposure conditions provided in the problem
light_level_lux = 50
# We assume a standard 8-hour day for an object "on display".
daily_exposure_hours = 8
days_per_year = 365

# --- Step 2: Calculate the annual light exposure ---
# Annual exposure = lux * hours_per_day * days_per_year
annual_exposure_lux_hours = light_level_lux * daily_exposure_hours * days_per_year

# --- Step 3: Calculate the time in years to reach the fading threshold ---
# Years to fade = Total JNF threshold / Annual exposure
time_to_fade_years = total_lux_hours_for_jnf / annual_exposure_lux_hours

# --- Step 4: Print the final equation and the result ---
print("The time in years for a just noticeable fade can be calculated as follows:")
print("Years to Fade = (Total Allowable Lux-Hours for JNF) / (Lux Level * Hours per Day * Days per Year)\n")
print("Using the standard values, the final equation is:")

# We format the numbers for clear presentation in the final equation.
# The result is rounded to two decimal places for readability.
print(f"{time_to_fade_years:.2f} years = {total_lux_hours_for_jnf} / ({light_level_lux} * {daily_exposure_hours} * {days_per_year})")
