import math

# --- Constants based on conservation science principles ---

# Total light exposure in lux-hours to cause a Just Noticeable Fade (JNF)
# for a highly sensitive material (ISO Bluewool 1) under damaging UV-rich light.
total_exposure_limit_lux_hours = 20000

# Given light intensity in lux.
light_intensity_lux = 50

# Assuming a standard museum display of 8 hours per day.
hours_per_day = 8

# Average number of days in a year to account for leap years.
days_per_year = 365.25

# --- Calculation ---

# Step 1: Calculate the total daily light exposure in lux-hours per day.
daily_exposure_lux_hours = light_intensity_lux * hours_per_day

# Step 2: Calculate the total number of days until a JNF occurs.
days_to_fade = total_exposure_limit_lux_hours / daily_exposure_lux_hours

# Step 3: Convert the number of days into years.
years_to_fade = days_to_fade / days_per_year

# --- Output ---

# Print the full equation showing all the numbers used in the calculation.
print(f"Time to Fade (years) = (Total Exposure Limit / (Light Intensity * Hours per Day)) / Days per Year")
print(f"Years = ({total_exposure_limit_lux_hours} / ({light_intensity_lux} * {hours_per_day})) / {days_per_year}")
print(f"Result: {years_to_fade:.3f} years")