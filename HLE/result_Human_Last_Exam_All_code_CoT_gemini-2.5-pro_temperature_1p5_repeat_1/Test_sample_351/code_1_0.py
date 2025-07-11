# Patient information
bsa_m2 = 0.95  # Body Surface Area in square meters

# Constants for calculation
standard_rate_ml_per_m2_per_day = 1500
hours_in_day = 24

# --- Calculation Steps ---

# 1. Calculate total daily maintenance fluid in mL
total_daily_fluid_ml = standard_rate_ml_per_m2_per_day * bsa_m2

# 2. Calculate the hourly rate in cc/hr (since 1 mL = 1 cc)
hourly_rate_cc_per_hr = total_daily_fluid_ml / hours_in_day

# --- Output the results ---

print("Calculating pediatric maintenance fluid needs using the Body Surface Area (BSA) method.")
print("\nPatient BSA: {} m^2".format(bsa_m2))

# Print the final equation with all the numbers
print("\nFinal Equation:")
print("({} mL/m²/day * {} m²) / {} hours/day = {:.2f} cc/hr".format(
    standard_rate_ml_per_m2_per_day,
    bsa_m2,
    hours_in_day,
    hourly_rate_cc_per_hr
))

print("\nThe required maintenance fluid rate is: {:.2f} cc/hr".format(hourly_rate_cc_per_hr))

# The final answer in the requested format
# Note: The result is rounded to two decimal places.
final_answer = "{:.2f}".format(hourly_rate_cc_per_hr)
print(f"\n<<<{final_answer}>>>")