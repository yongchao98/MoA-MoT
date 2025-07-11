# Patient and formula parameters
body_surface_area = 0.95  # in m^2
fluid_per_sq_meter_per_day = 1500  # in cc/m^2/day
hours_in_day = 24

# --- Calculations ---

# 1. Calculate total daily maintenance fluid needs
total_daily_fluid = fluid_per_sq_meter_per_day * body_surface_area

# 2. Calculate hourly maintenance fluid rate
hourly_fluid_rate = total_daily_fluid / hours_in_day

# --- Output the result ---

print("The formula for daily maintenance fluid needs based on Body Surface Area (BSA) is:")
print(f"Fluid Rate (cc/m^2/day) * BSA (m^2)\n")

print("Step 1: Calculate total daily fluid needs.")
print(f"Equation: {fluid_per_sq_meter_per_day} cc/m^2/day * {body_surface_area} m^2 = {total_daily_fluid:.1f} cc/day\n")

print("Step 2: Convert daily needs to hourly needs.")
print(f"Equation: {total_daily_fluid:.1f} cc/day / {hours_in_day} hours = {hourly_fluid_rate:.3f} cc/hr\n")

print(f"The calculated maintenance fluid rate for the patient is {hourly_fluid_rate:.3f} cc/hr.")