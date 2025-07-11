# Patient data
bsa = 0.95  # Body Surface Area in m^2
daily_fluid_rate_per_bsa = 1500 # Standard rate in mL/m^2/day

# Step 1: Explain the formula
print("The formula for calculating daily maintenance fluid needs based on Body Surface Area (BSA) is:")
print(f"Daily Fluid (mL) = {daily_fluid_rate_per_bsa} mL/m² * BSA (m²)")
print("-" * 30)

# Step 2: Calculate the total daily fluid requirement
total_daily_fluid = daily_fluid_rate_per_bsa * bsa
print("Calculating the total daily fluid volume:")
print(f"Daily Fluid = {daily_fluid_rate_per_bsa} mL/m² * {bsa} m² = {total_daily_fluid} mL/day")
print("-" * 30)

# Step 3: Calculate the hourly rate
hourly_rate = total_daily_fluid / 24
print("Calculating the hourly rate (cc/hr):")
print(f"Hourly Rate = {total_daily_fluid} mL / 24 hours = {hourly_rate:.1f} cc/hr")
print("-" * 30)

print(f"The final calculated maintenance fluid rate is {hourly_rate:.1f} cc/hr.")