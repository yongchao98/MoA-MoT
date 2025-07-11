import math

# --- Given Information ---
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500  # ml/day

# --- Step 1: Calculate Total Daily Fluid Requirement (Holliday-Segar Method) ---
if weight_kg <= 10:
    total_maintenance_fluid_day = weight_kg * 100
    calc_str_1 = f"{weight_kg} kg * 100 ml/kg = {total_maintenance_fluid_day} ml"
    hs_calculation_str = f"({calc_str_1})"
elif weight_kg <= 20:
    total_maintenance_fluid_day = (10 * 100) + ((weight_kg - 10) * 50)
    calc_str_1 = "10 kg * 100 ml/kg = 1000 ml"
    calc_str_2 = f"({weight_kg} - 10) kg * 50 ml/kg = {(weight_kg - 10) * 50} ml"
    hs_calculation_str = f"({calc_str_1}) + ({calc_str_2})"
else: # weight_kg > 20
    first_10_kg_fluid = 10 * 100
    next_10_kg_fluid = 10 * 50
    remaining_weight = weight_kg - 20
    remaining_fluid = remaining_weight * 20
    total_maintenance_fluid_day = first_10_kg_fluid + next_10_kg_fluid + remaining_fluid
    calc_str_1 = "10 kg * 100 ml/kg = 1000 ml"
    calc_str_2 = "10 kg * 50 ml/kg = 500 ml"
    calc_str_3 = f"({weight_kg} - 20) kg * 20 ml/kg = {remaining_fluid} ml"
    hs_calculation_str = f"({calc_str_1}) + ({calc_str_2}) + ({calc_str_3})"

print(f"Step 1: Calculate total daily fluid requirement for a {weight_kg} kg child.")
print(f"Holliday-Segar Calculation: {hs_calculation_str}")
print(f"Total Daily Fluid Requirement = {int(total_maintenance_fluid_day)} ml/day\n")


# --- Step 2: Calculate Daily Chemotherapy Fluid Volume ---
total_drug_dose_mg_day = drug_dose_per_m2_per_day * bsa_m2
chemo_fluid_volume_day = total_drug_dose_mg_day / admin_concentration

print("Step 2: Calculate the fluid volume from the chemotherapy drug per day.")
print(f"Daily Drug Dose (mg) = {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_drug_dose_mg_day} mg")
print(f"Daily Chemo Fluid Volume (ml) = {total_drug_dose_mg_day} mg / {admin_concentration} mg/ml = {int(chemo_fluid_volume_day)} ml\n")


# --- Step 3: Sum all non-maintenance fluid sources ---
total_other_fluids_day = chemo_fluid_volume_day + enteral_feeding_ml_day

print("Step 3: Sum all fluids from other sources.")
print(f"Total Other Fluids = {int(chemo_fluid_volume_day)} ml (chemo) + {enteral_feeding_ml_day} ml (milk) = {int(total_other_fluids_day)} ml/day\n")


# --- Step 4: Calculate Remaining Fluid Volume for Maintenance ---
iv_maintenance_fluid_day = total_maintenance_fluid_day - total_other_fluids_day

print("Step 4: Calculate the volume needed for IV maintenance fluid.")
print(f"IV Maintenance Fluid (ml/day) = {int(total_maintenance_fluid_day)} ml - {int(total_other_fluids_day)} ml = {int(iv_maintenance_fluid_day)} ml/day\n")


# --- Step 5: Calculate the Hourly Rate and Round ---
iv_rate_ml_hr = iv_maintenance_fluid_day / 24
rounded_iv_rate_ml_hr = round(iv_rate_ml_hr)

print("Step 5: Calculate the hourly rate for the maintenance fluid.")
print(f"IV Maintenance Rate (ml/hr) = {int(iv_maintenance_fluid_day)} ml / 24 hr = {iv_rate_ml_hr:.1f} ml/hr")
print(f"Rounded Rate = {rounded_iv_rate_ml_hr} ml/hr")

print("\n--- Final Answer ---")
print(f"The rate of the maintenance fluid should be set to {rounded_iv_rate_ml_hr} ml/hr.")