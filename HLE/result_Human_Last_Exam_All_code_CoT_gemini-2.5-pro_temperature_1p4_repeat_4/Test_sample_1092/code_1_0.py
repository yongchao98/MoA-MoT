import math

# Step 1: Define the standard fading dose for ISO Bluewool 1.
# This is the total exposure in lux-hours to cause a "Just Noticeable Fade" (JNF).
standard_fading_dose_lux_hours = 75000

# Step 2: Adjust the dose for the specified lighting conditions.
# "UV-rich light" can double the rate of fading. We represent this by
# halving the total exposure needed to cause a fade.
uv_factor = 2
adjusted_fading_dose = standard_fading_dose_lux_hours / uv_factor

# Step 3: Define the daily and annual exposure conditions.
lux_level = 50  # lux
hours_per_day = 8  # Assuming a standard 8-hour display day
days_per_year = 365 # Days in a year

# Calculate the total annual exposure in lux-hours per year.
annual_exposure = lux_level * hours_per_day * days_per_year

# Step 4: Calculate the number of years until the next JNF occurs.
years_to_fade = adjusted_fading_dose / annual_exposure

# Print the final calculation and result, showing each number in the equation.
print(f"To find the time to the next just noticeable fade, we perform the following calculation:")
print(f"Years to Fade = (Fading Dose for Bluewool 1 / UV Factor) / (Lux Level * Hours per Day * Days per Year)")
print(f"Years to Fade = ({standard_fading_dose_lux_hours} / {uv_factor}) / ({lux_level} * {hours_per_day} * {days_per_year})")
print(f"Years to Fade = {adjusted_fading_dose} / {annual_exposure}")
print(f"Result: {years_to_fade:.2f} years")
print(f"\nThis is approximately {math.ceil(years_to_fade * 12)} months.")
