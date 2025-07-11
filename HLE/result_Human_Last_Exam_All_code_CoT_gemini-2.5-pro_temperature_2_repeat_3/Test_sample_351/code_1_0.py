import math

def calculate_pediatric_burn_maintenance_fluids():
    """
    Calculates the maintenance fluid needs for a pediatric burn patient.

    This function combines the standard Holliday-Segar maintenance fluid calculation
    with the calculation for evaporative water loss (EWL) specific to burn patients.
    """
    # Patient data
    weight_kg = 25
    tbsa_percent = 45
    bsa_m2 = 0.95

    # --- Step 1: Calculate Standard Maintenance Fluids (Holliday-Segar Method) ---
    print("Step 1: Calculate baseline maintenance fluid needs using the Holliday-Segar method.")
    
    daily_maintenance_fluid = 0
    remaining_weight = weight_kg

    # For the first 10 kg
    if remaining_weight > 0:
        first_10_kg_fluid = min(remaining_weight, 10) * 100
        daily_maintenance_fluid += first_10_kg_fluid
        remaining_weight -= 10
        print(f"For the first 10 kg: 10 kg * 100 mL/kg = {int(first_10_kg_fluid)} mL")

    # For the next 10 kg (11-20 kg)
    if remaining_weight > 0:
        next_10_kg_fluid = min(remaining_weight, 10) * 50
        daily_maintenance_fluid += next_10_kg_fluid
        remaining_weight -= 10
        print(f"For the next 10 kg: 10 kg * 50 mL/kg = {int(next_10_kg_fluid)} mL")

    # For the remaining weight (>20 kg)
    if remaining_weight > 0:
        last_portion_fluid = remaining_weight * 20
        daily_maintenance_fluid += last_portion_fluid
        print(f"For the remaining {int(remaining_weight)} kg: {int(remaining_weight)} kg * 20 mL/kg = {int(last_portion_fluid)} mL")

    print(f"Total daily maintenance fluid: 1000 + 500 + 100 = {int(daily_maintenance_fluid)} mL/day")

    hourly_maintenance_rate = daily_maintenance_fluid / 24
    print(f"Hourly maintenance rate: {int(daily_maintenance_fluid)} mL / 24 hr = {hourly_maintenance_rate:.2f} cc/hr\n")


    # --- Step 2: Calculate Evaporative Water Loss (EWL) ---
    print("Step 2: Calculate additional fluid needs from evaporative water loss (EWL).")
    print("Using the formula: EWL (cc/hr) = (25 + %TBSA) * BSA")
    ewl_rate = (25 + tbsa_percent) * bsa_m2
    print(f"EWL = (25 + {tbsa_percent}) * {bsa_m2} = {ewl_rate:.2f} cc/hr\n")

    # --- Step 3: Calculate Total Maintenance Fluid Needs ---
    print("Step 3: Calculate total maintenance fluid needs by adding both rates.")
    total_fluid_rate = hourly_maintenance_rate + ewl_rate
    print("Total Rate (cc/hr) = Hourly Maintenance Rate + EWL Rate")
    print(f"Total Rate = {hourly_maintenance_rate:.2f} + {ewl_rate:.2f} = {total_fluid_rate:.2f} cc/hr\n")
    
    # Final rounded answer
    final_answer = round(total_fluid_rate, 1)
    print(f"The calculated total maintenance fluid need for the patient is approximately {final_answer} cc/hr.")
    
    return final_answer

# Execute the function and capture the final answer for submission format
final_result = calculate_pediatric_burn_maintenance_fluids()

print(f"\n<<<{final_result}>>>")