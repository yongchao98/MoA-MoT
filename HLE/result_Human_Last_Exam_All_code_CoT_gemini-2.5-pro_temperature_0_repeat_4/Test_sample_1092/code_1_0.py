# Plan:
# 1. Define the standard light dose (in lux-hours) for a Just Noticeable Fade (JNF) for an ISO Bluewool 1 object.
# 2. Adjust this dose for UV-rich light, which doubles the fading speed (halving the required dose).
# 3. Define the daily light exposure dose (in lux-hours).
# 4. Calculate the total annual light exposure dose.
# 5. Divide the adjusted JNF dose by the annual dose to find the time in years.
# 6. Print the full equation and the final result.

# Step 1: Standard JNF dose for Bluewool 1
standard_jnf_dose_bw1 = 50000  # in lux-hours

# Step 2: Adjust for UV-rich light (doubles the damage, so half the dose is needed)
uv_damage_factor = 2
effective_jnf_dose = standard_jnf_dose_bw1 / uv_damage_factor

# Step 3: Daily light exposure
daily_lux_hour_exposure = 50  # lux-hours per day

# Step 4: Calculate annual exposure
days_in_year = 365
annual_exposure = daily_lux_hour_exposure * days_in_year

# Step 5: Calculate years until the next just noticeable fade
years_to_fade = effective_jnf_dose / annual_exposure

# Step 6: Print the final equation and the result
# The final equation is: Time = Effective JNF Dose / (Daily Dose * Days per Year)
print(f"The time to the next just noticeable fade is calculated as follows:")
print(f"({standard_jnf_dose_bw1} lux-hours / {uv_damage_factor}) / ({daily_lux_hour_exposure} lux-hours/day * {days_in_year} days/year) = {years_to_fade:.2f} years")