import math

def calculate_maintenance_fluid_rate():
    """
    This function calculates the maintenance fluid rate for a child on chemotherapy
    based on the provided clinical data.
    """
    # Patient and drug information
    weight_kg = 22
    bsa_m2 = 0.8
    drug_dose_per_m2_per_day = 25  # mg/m^2/day
    admin_concentration_mg_per_ml = 1  # mg/ml
    enteral_feeding_ml_per_day = 500

    print("--- Step 1: Calculate the daily drug dose (mg) ---")
    daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
    print(f"Daily Dose = {drug_dose_per_m2_per_day} mg/m²/day * {bsa_m2} m² = {daily_drug_dose_mg} mg/day\n")

    print("--- Step 2: Calculate the daily fluid volume from the drug administration (ml) ---")
    daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration_mg_per_ml
    print(f"Drug Volume = {daily_drug_dose_mg} mg/day / {admin_concentration_mg_per_ml} mg/ml = {daily_drug_volume_ml} ml/day\n")

    print("--- Step 3: Calculate total daily maintenance fluid requirement (Holliday-Segar method) ---")
    if weight_kg <= 10:
        maintenance_fluid_per_day = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_fluid_per_day = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        fluid_first_10 = 10 * 100
        fluid_second_10 = 10 * 50
        fluid_remainder = (weight_kg - 20) * 20
        maintenance_fluid_per_day = fluid_first_10 + fluid_second_10 + fluid_remainder
    
    print(f"For a weight of {weight_kg} kg:")
    print(f"Requirement for first 10 kg: 10 * 100 = 1000 ml")
    print(f"Requirement for next 10 kg (11-20kg): 10 * 50 = 500 ml")
    print(f"Requirement for remaining {weight_kg - 20} kg: {weight_kg - 20} * 20 = {(weight_kg - 20) * 20} ml")
    print(f"Total Daily Requirement = 1000 + 500 + {(weight_kg - 20) * 20} = {maintenance_fluid_per_day} ml/day\n")

    print("--- Step 4: Calculate total fluid from other sources (drug + milk) ---")
    total_other_fluids_ml_per_day = daily_drug_volume_ml + enteral_feeding_ml_per_day
    print(f"Other Fluids = {daily_drug_volume_ml} ml (drug) + {enteral_feeding_ml_per_day} ml (milk) = {total_other_fluids_ml_per_day} ml/day\n")

    print("--- Step 5: Calculate the hourly rate of the maintenance fluid ---")
    remaining_maintenance_fluid_ml_per_day = maintenance_fluid_per_day - total_other_fluids_ml_per_day
    hourly_rate = remaining_maintenance_fluid_ml_per_day / 24
    rounded_hourly_rate = round(hourly_rate)

    print(f"The final equation for the hourly rate is:")
    print(f"Rate = round( (Total Requirement - Other Fluids) / 24 )")
    print(f"Rate = round( ({maintenance_fluid_per_day} - ({daily_drug_volume_ml} + {enteral_feeding_ml_per_day})) / 24 )")
    print(f"Rate = round( ({maintenance_fluid_per_day} - {total_other_fluids_ml_per_day}) / 24 )")
    print(f"Rate = round( {remaining_maintenance_fluid_ml_per_day} / 24 )")
    print(f"Rate = round( {hourly_rate:.2f} ) = {rounded_hourly_rate} ml/hr\n")
    
    print(f"The final calculated rate of the maintenance fluid is {rounded_hourly_rate} ml/hr.")
    return rounded_hourly_rate

# Run the calculation and store the final answer
final_answer = calculate_maintenance_fluid_rate()
# The final answer is wrapped in <<<>>>
print(f"\n<<<{final_answer}>>>")
