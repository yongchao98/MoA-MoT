import math

# Patient data
weight_kg = 25.0
tbsa_percent = 45.0
bsa_m2 = 0.95

# --- Step 1: Calculate Standard Maintenance Fluid (Holliday-Segar Method) ---
# Calculate daily maintenance fluid based on weight
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg
# 20 mL/kg for the remainder
first_10_kg_fluid = 10 * 100
second_10_kg_fluid = 10 * 50
remaining_weight = weight_kg - 20
remaining_fluid = remaining_weight * 20
daily_maintenance_ml = first_10_kg_fluid + second_10_kg_fluid + remaining_fluid

# Convert daily maintenance to an hourly rate
holliday_segar_rate_hr = daily_maintenance_ml / 24

print("--- Calculating Total Maintenance Fluid Needs ---")
print("\nThe calculation is composed of two parts: standard maintenance and evaporative loss.")

print("\n--- Part 1: Standard Maintenance Fluid (Holliday-Segar) ---")
print(f"For a patient weighing {weight_kg} kg:")
print(f"Calculation for daily fluid: (10kg * 100cc) + (10kg * 50cc) + ({remaining_weight}kg * 20cc) = {daily_maintenance_ml:.0f} cc/day")
print(f"Standard hourly maintenance rate = {daily_maintenance_ml:.0f} cc / 24 hr = {holliday_segar_rate_hr:.2f} cc/hr")

# --- Step 2: Calculate Evaporative Fluid Loss ---
# Use the formula: (25 + %TBSA) * BSA
evaporative_loss_rate_hr = (25 + tbsa_percent) * bsa_m2

print("\n--- Part 2: Additional Evaporative Loss from Burn ---")
print(f"For a {tbsa_percent}% TBSA burn and {bsa_m2} m^2 BSA:")
print(f"Calculation for hourly evaporative loss: (25 + {tbsa_percent}) * {bsa_m2} = {evaporative_loss_rate_hr:.2f} cc/hr")

# --- Step 3: Calculate Total Maintenance Fluid Rate ---
total_maintenance_rate_hr = holliday_segar_rate_hr + evaporative_loss_rate_hr

print("\n--- Part 3: Total Maintenance Fluid Rate ---")
print("Total Rate = Standard Hourly Rate + Hourly Evaporative Loss")
print(f"Final Equation: {holliday_segar_rate_hr:.2f} cc/hr + {evaporative_loss_rate_hr:.2f} cc/hr = {total_maintenance_rate_hr:.2f} cc/hr")

print(f"\nThe final calculated maintenance fluid rate is {total_maintenance_rate_hr:.1f} cc/hr.")
<<<133.2>>>