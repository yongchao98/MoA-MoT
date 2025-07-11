import math

# Step 1: Define the known values and assumptions.
# Standard light dose in lux-hours for a Just Noticeable Fade (JNF) for ISO Bluewool 1.
jnf_lux_hours_standard = 75000

# The light level the object is exposed to in lux.
lux_level = 50

# Step 2: Define factors based on conditions.
# Assumed damage factor for UV-rich light compared to the standard. A factor of 3 is a conservative estimate.
uv_damage_factor = 3

# Assumed hours of light exposure per day. A typical museum day is 8 hours.
hours_per_day = 8

# Days in a year, using 365.25 to account for leap years.
days_per_year = 365.25

# Step 3: Perform the calculations.
# Calculate the adjusted JNF dose, accounting for the damaging UV-rich light.
# The material is more sensitive, so it takes a smaller dose to fade.
adjusted_jnf_dose = jnf_lux_hours_standard / uv_damage_factor

# Calculate the total light dose the object receives per year.
annual_exposure = lux_level * hours_per_day * days_per_year

# Calculate the number of years until a JNF occurs.
years_to_fade = adjusted_jnf_dose / annual_exposure

# Step 4: Print the explanation and the result.
print("To calculate the years to the next just noticeable fade, we use the following formula:")
print("Years = (Standard JNF Dose / UV Damage Factor) / (Lux Level * Hours per Day * Days per Year)")
print("\nPlugging in the numbers:")
# The prompt requests each number in the final equation be output.
# Here we show the equation with all the defined values.
equation_string = (
    f"Years = ({jnf_lux_hours_standard} lux-hours / {uv_damage_factor}) / "
    f"({lux_level} lux * {hours_per_day} hours/day * {days_per_year} days/year)"
)
print(equation_string)

print(f"\nThis simplifies to:")
print(f"Years = {adjusted_jnf_dose:.0f} lux-hours / {annual_exposure:.0f} lux-hours/year")

print(f"\nResult:")
print(f"The next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")
<<<0.17>>>