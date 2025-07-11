# Define the known constants and assumptions
lux_hours_for_jnf_bw1 = 50000  # Total lux-hours for a JNF on Bluewool 1 material
lux_level = 50                  # Illuminance in lux
hours_per_day = 8               # Assumed hours of exposure per day
uv_damage_factor = 2            # Assumed damage factor for UV-rich light
days_in_year = 365.25           # Average days in a year to account for leap years

# Calculate the effective daily light exposure
effective_daily_exposure = lux_level * hours_per_day * uv_damage_factor

# Calculate the total number of days until a JNF occurs
days_to_fade = lux_hours_for_jnf_bw1 / effective_daily_exposure

# Convert the number of days to years
years_to_fade = days_to_fade / days_in_year

# Output the explanation and the final result
print("Calculation for years until next Just Noticeable Fade (JNF):")
print(f"Total Lux-Hours for JNF on Bluewool 1: {lux_hours_for_jnf_bw1}")
print(f"Daily exposure: {lux_level} lux * {hours_per_day} hours/day * {uv_damage_factor} (UV factor) = {effective_daily_exposure} effective lux-hours/day")
print(f"Years to fade = (Total Lux-Hours for JNF) / (Effective Daily Exposure * Days in Year)")
print(f"Years to fade = {lux_hours_for_jnf_bw1} / ({effective_daily_exposure} * {days_in_year})")
print(f"Result: {years_to_fade:.3f} years")
<<<0.171>>>