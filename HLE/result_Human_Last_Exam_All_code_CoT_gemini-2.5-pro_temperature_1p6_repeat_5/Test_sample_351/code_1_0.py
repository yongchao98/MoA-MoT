import math

# Patient variables
weight_kg = 25
tbsa_percent = 45

# --- Introduction of the Plan ---
print("To calculate the fluid needs for a pediatric burn patient, the standard of care is to combine")
print("the patient's normal daily maintenance fluid with the additional fluid required for burn resuscitation.")
print("The result is the total fluid volume needed over the first 24 hours, from which we can find the hourly rate.\n")

# --- Step 1: Calculate Standard Maintenance Fluid (Holliday-Segar Method) ---
print("Step 1: Calculate daily maintenance fluid using the Holliday-Segar formula.")
# The formula is:
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg
# 20 mL/kg for the remaining weight
part1 = 10 * 100
part2 = 10 * 50
part3 = (weight_kg - 20) * 20
maintenance_fluid_24hr = part1 + part2 + part3
print(f"For the first 10kg at 100 mL/kg: {part1} mL")
print(f"For the next 10kg at 50 mL/kg: {part2} mL")
print(f"For the remaining {weight_kg - 20}kg at 20 mL/kg: {part3} mL")
print(f"Total Daily Maintenance Fluid = {part1} + {part2} + {part3} = {maintenance_fluid_24hr} mL\n")

# --- Step 2: Calculate Burn Resuscitation Fluid (Parkland Formula) ---
print("Step 2: Calculate 24-hour burn resuscitation fluid using the Parkland formula (4 mL * kg * %TBSA).")
resuscitation_fluid_24hr = 4 * weight_kg * tbsa_percent
print(f"Burn Resuscitation Fluid = 4 mL * {weight_kg} kg * {tbsa_percent}% TBSA = {resuscitation_fluid_24hr} mL\n")

# --- Step 3: Calculate Total Fluid and Hourly Rate ---
print("Step 3: Combine fluids and calculate the average hourly rate.")
total_fluid_24hr = maintenance_fluid_24hr + resuscitation_fluid_24hr
hourly_rate = total_fluid_24hr / 24
print(f"Total 24-hour fluid = Maintenance ({maintenance_fluid_24hr} mL) + Resuscitation ({resuscitation_fluid_24hr} mL) = {total_fluid_24hr} mL")
print(f"Average hourly rate = Total Fluid ({total_fluid_24hr} mL) / 24 hours\n")

# --- Final Answer ---
print(f"The required maintenance fluid rate for this patient is {hourly_rate:.2f} cc/hr.")

<<<254.17>>>