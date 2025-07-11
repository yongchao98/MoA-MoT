import math

# --- Given Information ---
weight_kg = 22
bsa_m2 = 0.8
chemo_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# --- Step 1: Calculate Total Daily Maintenance Fluid (Holliday-Segar Method) ---
print("Step 1: Calculate Total Daily Maintenance Fluid Requirement (Holliday-Segar Method)")
if weight_kg <= 10:
    daily_fluid_req_ml = weight_kg * 100
    print(f"Fluid for {weight_kg} kg = {weight_kg} kg * 100 ml/kg = {daily_fluid_req_ml} ml/day")
else:
    fluid_for_first_10kg = 10 * 100
    print(f"Fluid for first 10 kg = 10 kg * 100 ml/kg = {fluid_for_first_10kg} ml")
    if weight_kg <= 20:
        fluid_for_next_kg = (weight_kg - 10) * 50
        print(f"Fluid for next {weight_kg - 10} kg = {weight_kg - 10} kg * 50 ml/kg = {fluid_for_next_kg} ml")
        daily_fluid_req_ml = fluid_for_first_10kg + fluid_for_next_kg
    else:
        fluid_for_next_10kg = 10 * 50
        fluid_for_remaining_kg = (weight_kg - 20) * 20
        print(f"Fluid for next 10 kg = 10 kg * 50 ml/kg = {fluid_for_next_10kg} ml")
        print(f"Fluid for remaining {weight_kg - 20} kg = {weight_kg - 20} kg * 20 ml/kg = {fluid_for_remaining_kg} ml")
        daily_fluid_req_ml = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg

print(f"Total Daily Fluid Requirement = {daily_fluid_req_ml} ml/day\n")


# --- Step 2: Calculate Daily Chemotherapy Volume ---
print("Step 2: Calculate Daily Chemotherapy Volume")
total_daily_dose_mg = bsa_m2 * chemo_dose_per_m2_per_day
print(f"Total daily drug dose = {bsa_m2} m^2 * {chemo_dose_per_m2_per_day} mg/m^2/day = {total_daily_dose_mg} mg/day")

chemo_volume_ml_day = total_daily_dose_mg / admin_concentration
print(f"Daily chemo volume = {total_daily_dose_mg} mg / {admin_concentration} mg/ml = {chemo_volume_ml_day} ml/day\n")


# --- Step 3: Sum Up All Non-Maintenance Fluid Volumes ---
print("Step 3: Sum Up All Non-Maintenance Fluid Volumes")
total_other_fluids_ml_day = enteral_feeding_ml_day + chemo_volume_ml_day
print(f"Total other fluids = {enteral_feeding_ml_day} ml (milk) + {chemo_volume_ml_day} ml (chemo) = {total_other_fluids_ml_day} ml/day\n")


# --- Step 4: Determine the Remaining Fluid Volume for Maintenance ---
print("Step 4: Determine the Remaining Fluid Volume for Maintenance")
maintenance_fluid_to_give_ml_day = daily_fluid_req_ml - total_other_fluids_ml_day
print(f"Remaining maintenance fluid = {daily_fluid_req_ml} ml - {total_other_fluids_ml_day} ml = {maintenance_fluid_to_give_ml_day} ml/day\n")


# --- Step 5: Calculate the Hourly Rate ---
print("Step 5: Calculate the Hourly Rate for Maintenance Fluid")
maintenance_rate_ml_hr = maintenance_fluid_to_give_ml_day / 24
print(f"Hourly rate = {maintenance_fluid_to_give_ml_day} ml / 24 hr = {maintenance_rate_ml_hr:.1f} ml/hr\n")


# --- Step 6: Round to the Nearest Whole Number ---
print("Step 6: Round to the Nearest Whole Number")
final_rate = round(maintenance_rate_ml_hr)
print(f"The rounded maintenance fluid rate is {final_rate} ml/hr.")

<<<43>>>