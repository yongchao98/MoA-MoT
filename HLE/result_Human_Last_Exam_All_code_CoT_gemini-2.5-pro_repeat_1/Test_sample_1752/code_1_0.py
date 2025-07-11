import math

def calculate_maintenance_fluid_rate():
    """
    This function calculates the maintenance fluid rate for a pediatric patient
    based on their weight, other fluid sources, and medication schedule.
    """
    # Patient and drug information
    weight_kg = 22
    bsa_m2 = 0.8
    drug_dose_per_m2 = 25  # mg/m²/day
    enteral_feeding_ml_day = 500
    drug_admin_concentration_mg_ml = 1.0

    # Step 1: Calculate total daily maintenance fluid requirement (Holliday-Segar method)
    print("Step 1: Calculate total daily maintenance fluid requirement (Holliday-Segar method).")
    print(f"Child's weight: {weight_kg} kg")
    
    if weight_kg <= 10:
        daily_fluid_needs_ml = weight_kg * 100
        print(f"Fluid for {weight_kg} kg: {weight_kg} kg * 100 ml/kg = {daily_fluid_needs_ml} ml")
    elif weight_kg <= 20:
        fluid_first_10 = 10 * 100
        fluid_next = (weight_kg - 10) * 50
        daily_fluid_needs_ml = fluid_first_10 + fluid_next
        print(f"Fluid for first 10 kg: 10 kg * 100 ml/kg = {fluid_first_10} ml")
        print(f"Fluid for next {weight_kg - 10} kg: {weight_kg - 10} kg * 50 ml/kg = {fluid_next} ml")
        print(f"Total Daily Fluid = {fluid_first_10} + {fluid_next} = {daily_fluid_needs_ml} ml")
    else: # weight_kg > 20
        fluid_first_10 = 10 * 100
        fluid_next_10 = 10 * 50
        remaining_weight = weight_kg - 20
        fluid_remaining = remaining_weight * 20
        daily_fluid_needs_ml = fluid_first_10 + fluid_next_10 + fluid_remaining
        print(f"Fluid for first 10 kg: 10 kg * 100 ml/kg = {fluid_first_10} ml")
        print(f"Fluid for next 10 kg: 10 kg * 50 ml/kg = {fluid_next_10} ml")
        print(f"Fluid for remaining {remaining_weight} kg: {remaining_weight} kg * 20 ml/kg = {fluid_remaining} ml")
        print(f"Total Daily Fluid = {fluid_first_10} + {fluid_next_10} + {fluid_remaining} = {daily_fluid_needs_ml} ml")
    print("-" * 30)

    # Step 2: Calculate the daily volume from the chemotherapy drug
    print("Step 2: Calculate the daily volume from the chemotherapy drug.")
    daily_drug_dose_mg = drug_dose_per_m2 * bsa_m2
    print(f"Daily drug dose = {drug_dose_per_m2} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg")
    
    daily_drug_volume_ml = daily_drug_dose_mg / drug_admin_concentration_mg_ml
    print(f"Daily drug volume = {daily_drug_dose_mg} mg / {drug_admin_concentration_mg_ml} mg/ml = {daily_drug_volume_ml} ml")
    print("-" * 30)

    # Step 3: Calculate the total fluid volume from other sources
    print("Step 3: Calculate the total fluid volume from other sources.")
    total_other_fluids_ml = daily_drug_volume_ml + enteral_feeding_ml_day
    print(f"Total from other sources = Drug volume + Enteral feeding")
    print(f"Total from other sources = {daily_drug_volume_ml} ml + {enteral_feeding_ml_day} ml = {total_other_fluids_ml} ml")
    print("-" * 30)

    # Step 4: Calculate the required volume for maintenance fluid
    print("Step 4: Calculate the required volume for maintenance fluid.")
    maintenance_fluid_volume_ml = daily_fluid_needs_ml - total_other_fluids_ml
    print(f"Maintenance fluid volume = Total Daily Fluid - Total from other sources")
    print(f"Maintenance fluid volume = {daily_fluid_needs_ml} ml - {total_other_fluids_ml} ml = {maintenance_fluid_volume_ml} ml")
    print("-" * 30)

    # Step 5: Calculate the infusion rate in ml/hr
    print("Step 5: Calculate the infusion rate in ml/hr.")
    infusion_rate_ml_hr = maintenance_fluid_volume_ml / 24
    print(f"Infusion rate = {maintenance_fluid_volume_ml} ml / 24 hours = {infusion_rate_ml_hr:.1f} ml/hr")
    print("-" * 30)

    # Step 6: Round to the nearest whole number
    print("Step 6: Round to the nearest whole number.")
    final_rate = round(infusion_rate_ml_hr)
    print(f"The final rate for the maintenance fluid, rounded to the nearest whole number, is {final_rate} ml/hr.")
    
    return final_rate

# Run the calculation and print the final answer in the required format
final_answer = calculate_maintenance_fluid_rate()
print(f"<<<{final_answer}>>>")
