# Patient-specific variables
bsa_m2 = 0.95  # Body Surface Area in square meters

# Standard constants for the BSA method
fluid_rate_per_m2_per_day = 1500  # mL/m^2/day
hours_per_day = 24

# --- Calculations ---

# 1. Calculate the total daily maintenance fluid volume
daily_fluid_needs_ml = fluid_rate_per_m2_per_day * bsa_m2

# 2. Calculate the hourly maintenance fluid rate
hourly_fluid_needs_ml_hr = daily_fluid_needs_ml / hours_per_day

# --- Output the result ---

print("Calculating pediatric maintenance fluid needs using the Body Surface Area (BSA) method.")
print("Formula: (Fluid Rate per m^2 * BSA) / Hours per Day\n")

# Print the calculation with the numbers used
print("Calculation:")
print(f"({fluid_rate_per_m2_per_day} mL/m^2/day * {bsa_m2} m^2) / {hours_per_day} hours = {hourly_fluid_needs_ml_hr:.1f} cc/hr")

print(f"\nThe required maintenance fluid rate is {hourly_fluid_needs_ml_hr:.1f} cc/hr.")