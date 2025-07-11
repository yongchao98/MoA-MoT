import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1.0  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate the total daily dose of the drug
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"Step 1: Calculate the daily drug dose.")
print(f"Daily Dose (mg) = Recommended Dose (mg/m^2) * BSA (m^2)")
print(f"Daily Dose = {drug_dose_per_m2_per_day} mg/m^2 * {bsa_m2} m^2 = {daily_drug_dose_mg} mg")
print("-" * 30)

# Step 2: Calculate the daily volume of the drug
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration
print(f"Step 2: Calculate the daily volume of the drug to be administered.")
print(f"Drug Volume (ml) = Daily Dose (mg) / Administration Concentration (mg/ml)")
print(f"Drug Volume = {daily_drug_dose_mg} mg / {admin_concentration} mg/ml = {daily_drug_volume_ml} ml")
print("-" * 30)

# Step 3: Calculate total daily maintenance fluid using Holliday-Segar method
maintenance_fluid_day = 0
if weight_kg > 20:
    fluid_first_10kg = 10 * 100
    fluid_next_10kg = 10 * 50
    fluid_remaining_kg = (weight_kg - 20) * 20
    maintenance_fluid_day = fluid_first_10kg + fluid_next_10kg + fluid_remaining_kg
    calculation_str = f"({10} kg * {100} ml/kg) + ({10} kg * {50} ml/kg) + ({(weight_kg - 20)} kg * {20} ml/kg)"
    result_str = f"{fluid_first_10kg} ml + {fluid_next_10kg} ml + {fluid_remaining_kg} ml"
elif weight_kg > 10:
    fluid_first_10kg = 10 * 100
    fluid_remaining_kg = (weight_kg - 10) * 50
    maintenance_fluid_day = fluid_first_10kg + fluid_remaining_kg
    calculation_str = f"({10} kg * {100} ml/kg) + ({(weight_kg - 10)} kg * {50} ml/kg)"
    result_str = f"{fluid_first_10kg} ml + {fluid_remaining_kg} ml"
else:
    maintenance_fluid_day = weight_kg * 100
    calculation_str = f"{weight_kg} kg * {100} ml/kg"
    result_str = f"{maintenance_fluid_day} ml"

print(f"Step 3: Calculate total daily maintenance fluid (Holliday-Segar).")
print(f"Calculation: {calculation_str}")
print(f"Daily Maintenance Fluid = {result_str} = {maintenance_fluid_day} ml/day")
print("-" * 30)

# Step 4: Calculate total fluid from other sources
total_other_fluids_ml = daily_drug_volume_ml + enteral_feeding_ml_day
print(f"Step 4: Calculate total fluid from other sources.")
print(f"Other Fluids (ml) = Drug Volume (ml) + Enteral Feeding (ml)")
print(f"Other Fluids = {daily_drug_volume_ml} ml + {enteral_feeding_ml_day} ml = {total_other_fluids_ml} ml")
print("-" * 30)

# Step 5: Calculate the remaining fluid volume to be given as maintenance fluid
remaining_fluid_ml = maintenance_fluid_day - total_other_fluids_ml
print(f"Step 5: Calculate the remaining maintenance fluid volume needed.")
print(f"Remaining Fluid (ml) = Total Maintenance (ml) - Other Fluids (ml)")
print(f"Remaining Fluid = {maintenance_fluid_day} ml - {total_other_fluids_ml} ml = {remaining_fluid_ml} ml")
print("-" * 30)

# Step 6: Calculate the rate of the maintenance fluid
maintenance_rate_ml_hr = remaining_fluid_ml / 24
rounded_rate = round(maintenance_rate_ml_hr)
print(f"Step 6: Calculate the hourly rate of the maintenance fluid.")
print(f"Rate (ml/hr) = Remaining Fluid (ml) / 24 hr")
print(f"Rate = {remaining_fluid_ml} ml / 24 hr = {maintenance_rate_ml_hr:.1f} ml/hr")
print(f"Rounded to the nearest whole number, the final rate is {rounded_rate} ml/hr.")

print(f"<<<{rounded_rate}>>>")