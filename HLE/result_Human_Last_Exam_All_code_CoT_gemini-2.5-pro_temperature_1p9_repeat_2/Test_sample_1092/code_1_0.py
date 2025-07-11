# Constants based on museum conservation standards
# For ISO Bluewool 1 under UV-rich light, one "Just Noticeable Fade" (JNF) occurs after approximately 50,000 lux-hours.
lux_hours_for_jnf = 50000

# Exposure conditions given in the problem
lux_level = 50  # lux
hours_per_day = 8  # Assuming a standard 8-hour exhibition day
days_per_year = 365.25  # Using 365.25 to account for leap years

# 1. Calculate the total exposure in lux-hours per year
total_lux_hours_per_year = lux_level * hours_per_day * days_per_year

# 2. Calculate the number of years until one JNF occurs
years_to_jnf = lux_hours_for_jnf / total_lux_hours_per_year

# 3. Print the full equation and the result
print(f"The time to the next just noticeable fade is calculated as follows:")
print(f"({lux_hours_for_jnf} lux-hours) / ({lux_level} lux * {hours_per_day} hours/day * {days_per_year} days/year) = {years_to_jnf:.2f} years")
print(f"\nAt this rate, the next just noticeable fade will occur in approximately {years_to_jnf:.2f} years.")
