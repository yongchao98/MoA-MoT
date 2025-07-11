import math

# --- Step 1: Define Constants ---
# Total light dose for one Just Noticeable Fade (JNF) for a Bluewool 1 rated material.
# This is a standard value from conservation science.
dose_for_jnf = 200000  # in lux-hours

# Lighting conditions from the user's query.
illuminance = 50  # in lux

# Assumption for "daily" exposure in a museum/display context.
hours_per_day = 8   # in hours
days_per_year = 365 # in days

# --- Step 2: Calculate Annual Dose ---
# This is the total lux-hours of exposure the object receives in one year.
annual_dose = illuminance * hours_per_day * days_per_year

# --- Step 3: Calculate Time to Fade ---
# The total time in years until one JNF is reached.
time_in_years = dose_for_jnf / annual_dose

# --- Output the results ---
# The user wants to see the numbers in the final equation.
print("To find the time to the next just noticeable fade (JNF) in years, we use the formula:")
print("Time (Years) = Total JNF Dose (lux-hours) / Annual Dose (lux-hours/year)")
print("\nWhere Annual Dose = Illuminance (lux) * Hours per Day * Days per Year")
print("\nUsing the provided and standard values:")
print(f"Time (Years) = {dose_for_jnf} / ({illuminance} * {hours_per_day} * {days_per_year})")
print(f"Time (Years) = {dose_for_jnf} / {annual_dose}")
print(f"\nResult: {time_in_years:.2f} years")
print("\nNote: This is an estimate. Fading will occur faster under UV-rich light than under UV-filtered light of the same intensity.")