# Patient Data
weight_kg = 25
tbsa_percent = 45
bsa_m2 = 0.95

# Convert TBSA percentage to a decimal
tbsa_decimal = tbsa_percent / 100

# 1. Calculate Standard Maintenance Fluids
# Formula: 1500 mL/m^2/day * BSA
maintenance_fluids_daily = 1500 * bsa_m2

# 2. Calculate Evaporative Fluid Loss
# Formula: 2000 mL/m^2/day * Burned Area
# Burned Area = BSA * TBSA%
burned_area = bsa_m2 * tbsa_decimal
evaporative_loss_daily = 2000 * burned_area

# 3. Calculate Total Daily Fluid Needs
total_fluids_daily = maintenance_fluids_daily + evaporative_loss_daily

# 4. Convert to Hourly Rate (cc/hr)
# 1 cc = 1 mL
fluid_rate_hr = total_fluids_daily / 24

# Print the step-by-step calculation
print("Step 1: Calculate Standard Maintenance Fluids (mL/day)")
print(f"    Equation: 1500 * BSA")
print(f"    Calculation: 1500 * {bsa_m2} = {maintenance_fluids_daily:.2f} mL/day\n")

print("Step 2: Calculate Evaporative Loss (mL/day)")
print(f"    Equation: 2000 * (BSA * TBSA%)")
print(f"    Calculation: 2000 * ({bsa_m2} * {tbsa_decimal}) = 2000 * {burned_area:.4f} = {evaporative_loss_daily:.2f} mL/day\n")

print("Step 3: Calculate Total Daily Fluid Needs (mL/day)")
print(f"    Equation: Maintenance + Evaporative Loss")
print(f"    Calculation: {maintenance_fluids_daily:.2f} + {evaporative_loss_daily:.2f} = {total_fluids_daily:.2f} mL/day\n")

print("Step 4: Convert to Hourly Rate (cc/hr)")
print(f"    Equation: Total mL/day / 24")
print(f"    Calculation: {total_fluids_daily:.2f} / 24 = {fluid_rate_hr:.2f} cc/hr\n")

print(f"Final Answer: The required maintenance fluid rate is {fluid_rate_hr:.0f} cc/hr.")
<<<95>>>