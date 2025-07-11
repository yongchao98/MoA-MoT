# 1. Define the known values based on conservation standards and the problem statement.

# Total allowed exposure in lux-hours for a Just Noticeable Fade (JNF) on a highly sensitive material (ISO Bluewool 1).
total_allowed_exposure = 50000

# Illuminance in lux.
illuminance_lux = 50

# Assumed daily exposure in hours for a museum on display.
hours_per_day = 8

# Number of days in a year.
days_per_year = 365

# 2. Calculate the total annual exposure in lux-hours per year.
annual_exposure = illuminance_lux * hours_per_day * days_per_year

# 3. Calculate the time in years for a Just Noticeable Fade to occur.
years_to_fade = total_allowed_exposure / annual_exposure

# 4. Print the final equation and the result.
print("Equation to calculate years to Just Noticeable Fade (JNF):")
print(f"Years = Total Allowed Exposure / (Illuminance * Hours per Day * Days per Year)")
print(f"Years = {total_allowed_exposure} / ({illuminance_lux} * {hours_per_day} * {days_per_year})")
print(f"\nTime to next just noticeable fade: {years_to_fade:.2f} years.")
