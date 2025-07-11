# Patient's Body Surface Area in square meters
bsa = 0.95

# Standard fluid rate per square meter per day (in mL)
daily_rate_per_m2 = 1500

# 1. Calculate the total daily fluid requirement
total_daily_fluid = daily_rate_per_m2 * bsa

# 2. Convert the daily requirement to an hourly rate
hours_in_day = 24
hourly_rate = total_daily_fluid / hours_in_day

# 3. Print the final equation and the answer
# Note: cc (cubic centimeter) is equivalent to mL (milliliter)
print("Calculation for Maintenance Fluid Rate (cc/hr):")
print(f"({daily_rate_per_m2} mL/m^2/day * {bsa} m^2) / {hours_in_day} hours = {hourly_rate:.2f} cc/hr")
