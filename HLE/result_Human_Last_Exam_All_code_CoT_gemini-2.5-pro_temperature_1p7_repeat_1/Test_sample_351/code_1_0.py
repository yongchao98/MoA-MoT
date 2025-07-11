# Patient parameters
tbsa_percent = 45  # %
bsa_m2 = 0.95      # m^2

# Constants from the Galveston formula
fluid_rate_burn_loss = 5000  # mL/m^2 burned
fluid_rate_maintenance = 2000  # mL/m^2 total BSA
hours_in_day = 24

# --- Step-by-step Calculation ---

# 1. Calculate the surface area that is burned in m^2
bsa_burned_m2 = bsa_m2 * (tbsa_percent / 100.0)

# 2. Calculate the fluid volume for burn loss
fluid_for_burn_loss_daily = fluid_rate_burn_loss * bsa_burned_m2

# 3. Calculate the maintenance fluid volume
fluid_for_maintenance_daily = fluid_rate_maintenance * bsa_m2

# 4. Calculate the total daily fluid requirement
total_daily_fluid = fluid_for_burn_loss_daily + fluid_for_maintenance_daily

# 5. Calculate the hourly rate in cc/hr (1 mL = 1 cc)
hourly_rate_cc = total_daily_fluid / hours_in_day

# --- Output the results ---

print("Calculating pediatric burn fluid needs using the Galveston Formula.")
print(f"Formula: (({fluid_rate_burn_loss} mL/m^2 * Burned BSA) + ({fluid_rate_maintenance} mL/m^2 * Total BSA)) / {hours_in_day} hours")
print("-" * 30)
print(f"Given:")
print(f"  Total Body Surface Area (BSA): {bsa_m2} m^2")
print(f"  Burn Percentage (TBSA): {tbsa_percent}%")
print("-" * 30)

print("Step 1: Calculate Fluid for Burn Loss")
print(f"  Burned BSA = {tbsa_percent / 100.0} * {bsa_m2} m^2 = {bsa_burned_m2:.4f} m^2")
print(f"  Fluid for Burn Loss = {fluid_rate_burn_loss} mL/m^2 * {bsa_burned_m2:.4f} m^2 = {fluid_for_burn_loss_daily:.1f} mL/day")
print("\nStep 2: Calculate Fluid for Maintenance")
print(f"  Fluid for Maintenance = {fluid_rate_maintenance} mL/m^2 * {bsa_m2} m^2 = {fluid_for_maintenance_daily:.1f} mL/day")
print("\nStep 3: Calculate Total Hourly Rate")
print(f"  Total Daily Fluid = {fluid_for_burn_loss_daily:.1f} mL + {fluid_for_maintenance_daily:.1f} mL = {total_daily_fluid:.1f} mL/day")
print("-" * 30)
print("Final Calculation:")
print(f"({fluid_rate_burn_loss} * ({tbsa_percent} / 100) * {bsa_m2} + {fluid_rate_maintenance} * {bsa_m2}) / {hours_in_day} = {hourly_rate_cc:.1f} cc/hr")
print("-" * 30)
print(f"The required maintenance fluid rate is {hourly_rate_cc:.1f} cc/hr.")
