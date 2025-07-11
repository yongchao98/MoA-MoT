import math

# Step 1 & 2: Define constants and the UV-adjusted light damage budget.
# ISO Bluewool 1 corresponds to ~50,000 lux-hours for one Just Noticeable Fade (JNF).
jnf_limit_bw1 = 50000  # lux-hours

# UV-rich light can double the rate of damage, so we divide the limit by a factor of 2.
uv_damage_factor = 2
adjusted_jnf_limit = jnf_limit_bw1 / uv_damage_factor

# Step 3: Calculate the annual light exposure.
# Light intensity is 50 lux.
lux = 50
# We assume a standard museum day of 8 hours.
hours_per_day = 8
days_per_year = 365

annual_exposure = lux * hours_per_day * days_per_year

# Step 4: Calculate the number of years to reach one JNF.
years_to_fade = adjusted_jnf_limit / annual_exposure

# Output the equation with the numbers used in the calculation.
print(f"Equation: Years to Fade = (JNF Limit for BW1 / UV Damage Factor) / (Lux Level * Hours per Day * Days per Year)")
print(f"Calculation: Years to Fade = ({jnf_limit_bw1} / {uv_damage_factor}) / ({lux} * {hours_per_day} * {days_per_year})")
print(f"Result: {years_to_fade:.3f} years")
print("\nThis means a just noticeable fade will occur in approximately 0.171 years, which is about 2 months.")

# Final answer in the required format
final_answer = round(years_to_fade, 3)