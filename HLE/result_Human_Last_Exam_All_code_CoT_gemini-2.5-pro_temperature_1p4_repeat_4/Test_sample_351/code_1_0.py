import math

# Patient and formula constants
bsa_m2 = 0.95
fluid_per_bsa_per_day = 1500
hours_in_day = 24

# Step 1: Calculate total daily maintenance fluid.
# Formula: Daily Fluid (mL) = 1500 mL/m^2/day * BSA (m^2)
total_daily_fluid = fluid_per_bsa_per_day * bsa_m2

# Step 2: Calculate the hourly rate.
# Formula: Hourly Rate (mL/hr) = Total Daily Fluid / 24 hours
hourly_fluid_rate = total_daily_fluid / hours_in_day

# Output the results and the equation
print("Calculating pediatric maintenance fluid needs using the Body Surface Area (BSA) method.")
print(f"The formula is: (Fluid Rate per BSA * Patient's BSA) / Hours in a Day\n")

print(f"The final calculation is based on the following equation:")
# The user requested each number in the final equation be output
print(f"({fluid_per_bsa_per_day} mL/m²/day * {bsa_m2} m²) / {hours_in_day} hours = {hourly_fluid_rate:.2f} cc/hr")

print(f"\nThe required maintenance fluid rate for the patient is {hourly_fluid_rate:.2f} cc/hr.")