import math

# Patient Information
weight_kg = 25
tbsa_percent = 45
bsa_m2 = 0.95

# Constants for calculation
standard_rate_per_bsa_per_day = 1500  # mL/m^2/day
hours_in_day = 24

# This script calculates maintenance fluid needs using the Body Surface Area (BSA) method.
# This is a common method for pediatric patients and is suggested by the provision of BSA in the problem.

# Step 1: Calculate total daily maintenance fluid volume.
daily_fluid_volume = standard_rate_per_bsa_per_day * bsa_m2

# Step 2: Convert the daily volume to an hourly rate.
hourly_fluid_rate = daily_fluid_volume / hours_in_day

# Step 3: Print the calculation process and the final result.
print("Calculating pediatric maintenance fluid needs using the BSA method.")
print(f"Patient's BSA: {bsa_m2} m^2")
print("-" * 30)

print("1. Calculate total daily fluid volume (mL/day):")
print(f"   Formula: Standard Rate (mL/m^2/day) * BSA (m^2)")
print(f"   Calculation: {standard_rate_per_bsa_per_day} * {bsa_m2} = {daily_fluid_volume:.1f} mL/day\n")

print("2. Convert daily volume to hourly rate (cc/hr):")
print(f"   Formula: Daily Volume (mL/day) / 24 hours")
print(f"   Calculation: {daily_fluid_volume:.1f} / {hours_in_day} = {hourly_fluid_rate:.1f} cc/hr\n")

print(f"The final calculated maintenance fluid rate is approximately {hourly_fluid_rate:.1f} cc/hr.")

# The final answer is the numerical value of the hourly fluid rate.
final_answer = round(hourly_fluid_rate, 1)
# print(f"<<<{final_answer}>>>") # This would be for internal processing.