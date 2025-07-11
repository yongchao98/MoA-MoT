# Define the known constants and variables
# For ISO Bluewool 1, a Just Noticeable Fade (JNF) occurs at ~50,000 lux-hours.
total_lux_hours_for_jnf = 50000

# Light intensity in lux
lux_level = 50

# Assumed hours of exposure per day in a museum setting
hours_per_day = 8

# Days in a year, accounting for leap years
days_per_year = 365.25

# Step 1: Calculate the total annual light exposure in lux-hours
annual_exposure = lux_level * hours_per_day * days_per_year

# Step 2: Calculate the number of years until a Just Noticeable Fade occurs
years_to_fade = total_lux_hours_for_jnf / annual_exposure

# Print the final equation with all the numbers and the result
print(f"Time to Just Noticeable Fade Calculation:")
print(f"({total_lux_hours_for_jnf} lux-hours) / ({lux_level} lux * {hours_per_day} hours/day * {days_per_year} days/year) = {years_to_fade:.2f} years")