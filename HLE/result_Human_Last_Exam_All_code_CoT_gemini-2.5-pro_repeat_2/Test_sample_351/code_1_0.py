import math

# Patient and formula variables
weight_kg = 25
tbsa_percent = 45
bsa_m2 = 0.95
maintenance_factor_per_day = 1500 # mL per m^2 per day
hours_in_day = 24

# --- Maintenance Fluid Calculation (BSA Method) ---
# This is the primary request of the user.

# 1. Calculate total daily maintenance fluid
daily_maintenance_fluid = maintenance_factor_per_day * bsa_m2

# 2. Calculate hourly maintenance fluid rate
hourly_maintenance_rate = daily_maintenance_fluid / hours_in_day

# --- Output the result ---
print("Calculating pediatric maintenance fluid needs based on Body Surface Area (BSA).")
print(f"The formula is: (Maintenance Factor * BSA) / Hours in a Day")
print("\nStep 1: Calculate total daily fluid volume.")
print(f"Equation: {maintenance_factor_per_day} mL/m²/day * {bsa_m2} m² = {daily_maintenance_fluid:.2f} mL/day")

print("\nStep 2: Convert daily volume to an hourly rate.")
print(f"Equation: {daily_maintenance_fluid:.2f} mL/day / {hours_in_day} hours/day = {hourly_maintenance_rate:.2f} cc/hr")

print("\n---------------------------------------------------------")
print("Final Calculation:")
# The user requested to see the numbers in the final equation.
print(f"({maintenance_factor_per_day} * {bsa_m2}) / {hours_in_day} = {hourly_maintenance_rate:.2f} cc/hr")
print("---------------------------------------------------------")

# Note: The total fluid for a burn patient would also include resuscitation fluids (e.g., via Parkland formula),
# but this calculation correctly addresses the specific request for *maintenance fluid needs*.
final_answer = round(hourly_maintenance_rate, 2)