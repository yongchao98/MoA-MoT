# Step 1: Define the constants and variables based on the problem and conservation science standards.

# The total light exposure (in lux-hours) to cause a "Just Noticeable Fade" (JNF) for ISO Bluewool #1.
total_lux_hours_for_jnf = 20000

# The illuminance of the light source in lux.
lux_level = 50

# The number of hours the object is exposed to light per day. We assume a standard 8-hour day.
hours_per_day = 8

# The acceleration factor for UV-rich light compared to standard low-UV museum light. We assume a factor of 3.
uv_factor = 3

# The number of days in a year, accounting for leap years on average.
days_per_year = 365.25

# Step 2: Calculate the total number of years until a JNF occurs.

# Calculate the effective lux-hours of exposure per day.
effective_daily_lux_hours = lux_level * hours_per_day * uv_factor

# Calculate the total number of days to reach JNF.
total_days_to_jnf = total_lux_hours_for_jnf / effective_daily_lux_hours

# Convert the total number of days into years.
years_to_jnf = total_days_to_jnf / days_per_year

# Step 3: Print the final equation and the result.
# The equation shows how the total lux-hours for fading is divided by the annual effective light exposure.
print("Calculation for years until Just Noticeable Fade (JNF):")
print(f"({total_lux_hours_for_jnf} lux-hours for JNF) / (({lux_level} lux * {hours_per_day} hours/day * {uv_factor} UV-factor) * {days_per_year} days/year)")
print(f"= {years_to_jnf:.4f} years")
