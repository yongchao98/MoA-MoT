import math

# Step 1: Define constants and assumptions.
# The cumulative exposure to cause one Just Noticeable Fade (JNF) for ISO Bluewool 1.
jnf_lux_hours_bw1 = 50000
# The light level the object is exposed to.
lux_level = 50
# Assumption: 'daily' exposure means lights are on for 8 hours.
hours_per_day = 8
# Assumption: 'UV-rich' light is twice as damaging as standard UV-filtered light.
uv_damage_factor = 2
# Average number of days in a year to account for leap years.
days_per_year = 365.25

# Step 2: Calculate the effective daily exposure.
# The effective lux level considering the UV damage.
effective_lux = lux_level * uv_damage_factor
# The total effective exposure per day.
daily_effective_lux_hours = effective_lux * hours_per_day

# Step 3: Calculate the time until one JNF occurs.
# Total days until a JNF.
days_to_fade = jnf_lux_hours_bw1 / daily_effective_lux_hours
# Total years until a JNF.
years_to_fade = days_to_fade / days_per_year

# Step 4: Print the plan, the equation, and the result.
print("Calculating the years until the next Just Noticeable Fade (JNF).")
print("-" * 60)
print(f"Threshold for ISO Bluewool 1 (JNF): {jnf_lux_hours_bw1} lux hours")
print(f"Light level: {lux_level} lux")
print(f"Assumed daily exposure: {hours_per_day} hours")
print(f"Assumed UV damage factor: {uv_damage_factor}")
print("-" * 60)

print("The final equation is:")
print("Years = (Total JNF Lux Hours) / (Lux Level * UV Factor * Hours per Day * Days per Year)")
print(f"Years = {jnf_lux_hours_bw1} / ({lux_level} * {uv_damage_factor} * {hours_per_day} * {days_per_year})")
print(f"\nResult: The next just noticeable fade will occur in approximately {years_to_fade:.3f} years.")
