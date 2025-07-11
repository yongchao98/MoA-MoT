import math

def calculate_fluid_rate():
    """
    Calculates the maintenance fluid rate for a pediatric patient.
    """
    # Patient and prescription data
    weight_kg = 22
    bsa_m2 = 0.8
    drug_dose_per_m2_per_day = 25  # mg/m²/day
    drug_admin_concentration = 1  # mg/ml
    enteral_feeding_ml_per_day = 500

    # Step 1: Calculate Total Daily Maintenance Fluid (Holliday-Segar Method)
    print("--- Step 1: Calculate Total Daily Maintenance Fluid (Holliday-Segar) ---")
    fluid_first_10kg = 10 * 100
    fluid_next_10kg = 10 * 50
    fluid_remaining_kg = (weight_kg - 20) * 20
    daily_maintenance_fluid = fluid_first_10kg + fluid_next_10kg + fluid_remaining_kg
    
    print(f"For the first 10 kg: 10 kg * 100 ml/kg = {fluid_first_10kg} ml")
    print(f"For the next 10 kg: 10 kg * 50 ml/kg = {fluid_next_10kg} ml")
    print(f"For the remaining {weight_kg - 20} kg: ({weight_kg} - 20) kg * 20 ml/kg = {fluid_remaining_kg} ml")
    print(f"Total Daily Maintenance Fluid = {fluid_first_10kg} + {fluid_next_10kg} + {fluid_remaining_kg} = {daily_maintenance_fluid} ml/day\n")

    # Step 2: Calculate Daily Volume of Other Fluids
    print("--- Step 2: Calculate Daily Volume of Other Fluids ---")
    # Drug volume
    total_daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
    daily_drug_volume_ml = total_daily_drug_dose_mg / drug_admin_concentration
    print(f"Daily drug volume = ({drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m²) / {drug_admin_concentration} mg/ml = {daily_drug_volume_ml} ml")

    # Total other fluids
    total_other_fluids = daily_drug_volume_ml + enteral_feeding_ml_per_day
    print(f"Total volume from other sources = {daily_drug_volume_ml} ml (drug) + {enteral_feeding_ml_per_day} ml (feeding) = {total_other_fluids} ml/day\n")

    # Step 3: Calculate Remaining Volume for IV Maintenance Fluid
    print("--- Step 3: Calculate Remaining Volume for IV Maintenance Fluid ---")
    remaining_iv_fluid = daily_maintenance_fluid - total_other_fluids
    print(f"Remaining IV fluid volume = {daily_maintenance_fluid} ml - {total_other_fluids} ml = {remaining_iv_fluid} ml/day\n")
    
    # Step 4: Calculate and Finalize the Hourly Rate
    print("--- Step 4: Calculate and Finalize the Hourly Rate ---")
    hourly_rate = remaining_iv_fluid / 24
    final_rate = round(hourly_rate)
    
    print("Final Calculation:")
    print(f"round( ( (10 * 100 + 10 * 50 + (22 - 20) * 20) - ( (25 * 0.8) / 1 + 500 ) ) / 24 )")
    print(f"= round( ( ({fluid_first_10kg} + {fluid_next_10kg} + {fluid_remaining_kg}) - ({daily_drug_volume_ml} + {enteral_feeding_ml_per_day}) ) / 24 )")
    print(f"= round( ( {daily_maintenance_fluid} - {total_other_fluids} ) / 24 )")
    print(f"= round( {remaining_iv_fluid} / 24 )")
    print(f"= round( {hourly_rate:.2f} )")
    print(f"= {final_rate} ml/hr")
    
    return final_rate

if __name__ == '__main__':
    rate = calculate_fluid_rate()
    print(f"\nThe final maintenance fluid rate is: {rate} ml/hr.")
    print(f"<<<{rate}>>>")
