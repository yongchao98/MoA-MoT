import math

# Patient data
weight_kg = 25
tbsa_percent = 45
bsa_m2 = 0.95

# Step 1: Calculate daily maintenance fluid (Holliday-Segar method)
# For weight > 20kg
maint_fluid_part1 = 10 * 100
maint_fluid_part2 = 10 * 50
maint_fluid_part3 = (weight_kg - 20) * 20
total_daily_maint_fluid = maint_fluid_part1 + maint_fluid_part2 + maint_fluid_part3
maint_rate_hr = total_daily_maint_fluid / 24

# Step 2: Calculate hourly evaporative water loss (EWL)
ewl_rate_hr = (25 + tbsa_percent) * bsa_m2

# Step 3: Calculate total hourly fluid rate
total_hourly_rate = maint_rate_hr + ewl_rate_hr

# Displaying the calculation process and result
print("The calculation for total maintenance fluid in a pediatric burn patient involves adding baseline needs to evaporative loss.")
print("\n--- Calculation Breakdown ---")

print(f"\n1. Baseline Maintenance Fluid (Holliday-Segar):")
print(f"   Daily Total = ({10} kg * {100} mL/kg) + ({10} kg * {50} mL/kg) + ({(weight_kg - 20)} kg * {20} mL/kg) = {total_daily_maint_fluid} mL/day")
print(f"   Hourly Rate = {total_daily_maint_fluid} mL / 24 hr = {maint_rate_hr:.1f} cc/hr")


print(f"\n2. Evaporative Water Loss (EWL):")
print(f"   Hourly Rate = (25 + {tbsa_percent} % TBSA) * {bsa_m2} m^2 BSA = {ewl_rate_hr:.1f} cc/hr")

print("\n--- Final Calculation ---")
print("Total Fluid Rate = Baseline Rate + EWL Rate")
print(f"Total Fluid Rate (cc/hr) = {maint_rate_hr:.1f} + {ewl_rate_hr:.1f} = {total_hourly_rate:.1f}")

# Final result rounded to one decimal place
final_answer = round(total_hourly_rate, 1)
print(f"\nThe required maintenance fluid rate is {final_answer} cc/hr.")
