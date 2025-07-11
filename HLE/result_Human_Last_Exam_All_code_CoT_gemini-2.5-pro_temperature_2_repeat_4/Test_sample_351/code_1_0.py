import math

# Patient's given data
bsa = 0.95  # Body Surface Area in m^2
daily_rate_per_bsa = 1500 # Standard maintenance rate in mL/m^2/day

# Step 1: Calculate total daily maintenance fluid
daily_maintenance_fluid = daily_rate_per_bsa * bsa

# Step 2: Convert daily fluid requirement to an hourly rate
hourly_maintenance_rate = daily_maintenance_fluid / 24

# Print the final equation and the result
print("Calculation of Hourly Maintenance Fluid Rate:")
print(f"Daily Fluid Requirement (mL) = {daily_rate_per_bsa} mL/m^2/day * {bsa} m^2 = {daily_maintenance_fluid:.2f} mL/day")
print(f"Hourly Rate (cc/hr) = {daily_maintenance_fluid:.2f} mL/day / 24 hours = {hourly_maintenance_rate:.1f} cc/hr")

final_answer = round(hourly_maintenance_rate, 1)
# The final answer will be returned in the format <<<answer content>>> below
# print(f"<<<{final_answer}>>>")