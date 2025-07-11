import math

# Step 1: Define the constants and assumptions based on conservation science standards.

# The total cumulative light exposure that causes one "Just Noticeable Fade" (JNF)
# for a highly sensitive material rated as ISO Bluewool 1.
total_lux_hours_for_jnf = 1500000

# The intensity of the light source, as given in the problem.
lux_level = 50

# A damage factor for UV-rich light, which accelerates fading.
# It's estimated to be about 3 times more damaging than standard incandescent light.
uv_damage_factor = 3

# Assume a standard display schedule for a museum object.
hours_per_day = 8
days_per_year = 365

# Step 2: Calculate the total damaging exposure the object receives per year.
annual_damaging_exposure = lux_level * uv_damage_factor * hours_per_day * days_per_year

# Step 3: Calculate the number of years until one JNF occurs by dividing the
# JNF threshold by the annual damaging exposure.
years_to_fade = total_lux_hours_for_jnf / annual_damaging_exposure

# Step 4: Print the final equation with all the numerical values and the result.
print("This script calculates the time in years for a just noticeable fade to occur.")
print("---")
print("The final calculation is based on the following equation:")
print("Years to Fade = Total Lux-Hours for JNF / (Lux Level * UV Factor * Hours per Day * Days per Year)")
print("\nPlugging in the numbers:")
# The format string uses {:,} to add thousands separators for readability.
print(f"{total_lux_hours_for_jnf:,} / ({lux_level} * {uv_damage_factor} * {hours_per_day} * {days_per_year}) = {years_to_fade:.2f} years")
print("---")
print(f"It will take approximately {years_to_fade:.2f} years for the next just noticeable fade to occur.")
