# Patient's given data
# Body Surface Area in square meters (m^2)
bsa = 0.95
# Standard maintenance fluid rate in mL per square meter per day
maintenance_rate_per_m2_per_day = 1500
# Number of hours in a day
hours_in_day = 24

# --- Calculation ---

# 1. Calculate the total fluid needed over 24 hours
total_daily_fluid = bsa * maintenance_rate_per_m2_per_day

# 2. Calculate the hourly fluid rate in cc/hr (mL/hr)
hourly_rate = total_daily_fluid / hours_in_day

# --- Output ---
print("Calculation of Pediatric Maintenance Fluids using the BSA Method")
print("---------------------------------------------------------------")
print(f"Formula: (BSA * Fluid Rate per m^2) / Hours in a day")
print(f"The calculation is: ({bsa} m^2 * {maintenance_rate_per_m2_per_day} mL/m^2/day) / {hours_in_day} hours")
print(f"Result: The required maintenance fluid rate is {hourly_rate:.2f} cc/hr.")