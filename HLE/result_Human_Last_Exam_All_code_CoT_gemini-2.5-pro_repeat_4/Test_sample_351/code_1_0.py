import math

def calculate_pediatric_burn_fluid():
    """
    Calculates the total fluid needs for a pediatric burn patient using the
    Modified Parkland Formula (Parkland + Holliday-Segar Maintenance).
    """
    # Patient data
    weight_kg = 25
    tbsa_percent = 45

    # --- Step 1: Calculate 24-hour maintenance fluid (Holliday-Segar formula) ---
    # 100 mL/kg for the first 10 kg
    # 50 mL/kg for the next 10 kg (up to 20 kg)
    # 20 mL/kg for the remaining weight
    
    # These variables store the fluid for each weight bracket
    fluid_part1 = 0
    fluid_part2 = 0
    fluid_part3 = 0
    
    # These variables store the weight in each bracket for printing
    weight_part1 = 0
    weight_part2 = 0
    weight_part3 = 0

    if weight_kg > 20:
        weight_part1 = 10
        weight_part2 = 10
        weight_part3 = weight_kg - 20
        fluid_part1 = weight_part1 * 100
        fluid_part2 = weight_part2 * 50
        fluid_part3 = weight_part3 * 20
    elif weight_kg > 10:
        weight_part1 = 10
        weight_part2 = weight_kg - 10
        fluid_part1 = weight_part1 * 100
        fluid_part2 = weight_part2 * 50
    else:
        weight_part1 = weight_kg
        fluid_part1 = weight_part1 * 100

    maintenance_fluid_24hr = fluid_part1 + fluid_part2 + fluid_part3

    # --- Step 2: Calculate 24-hour burn resuscitation fluid (Parkland formula) ---
    # Formula: 4 mL * weight (kg) * % TBSA
    burn_fluid_24hr = 4 * weight_kg * tbsa_percent

    # --- Step 3: Calculate total fluid and average hourly rate ---
    total_fluid_24hr = maintenance_fluid_24hr + burn_fluid_24hr
    average_hourly_rate = total_fluid_24hr / 24

    # --- Print the detailed calculation ---
    print("Calculation of Total Fluid Requirement for a Pediatric Burn Patient")
    print("-" * 60)
    print(f"Patient Data: Weight = {weight_kg} kg, TBSA Burn = {tbsa_percent}%")
    print("-" * 60)

    print("Step 1: Calculate Daily Maintenance Fluid (Holliday-Segar Formula)")
    print(f"For first {weight_part1} kg: {weight_part1} kg * 100 mL/kg = {fluid_part1} mL")
    if weight_part2 > 0:
        print(f"For next {weight_part2} kg:  {weight_part2} kg * 50 mL/kg  = {fluid_part2} mL")
    if weight_part3 > 0:
        print(f"For remaining {weight_part3} kg: {weight_part3} kg * 20 mL/kg  = {fluid_part3} mL")
    print(f"Total Maintenance Fluid = {fluid_part1} + {fluid_part2} + {fluid_part3} = {maintenance_fluid_24hr} mL/day")
    print("-" * 60)

    print("Step 2: Calculate Burn Resuscitation Fluid (Parkland Formula)")
    print(f"Formula: 4 mL * Weight (kg) * %TBSA")
    print(f"Burn Fluid = 4 * {weight_kg} * {tbsa_percent} = {burn_fluid_24hr} mL/day")
    print("-" * 60)

    print("Step 3: Calculate Total Fluid and Average Hourly Rate")
    print(f"Total Fluid (24hr) = Maintenance Fluid + Burn Fluid")
    print(f"Total Fluid (24hr) = {maintenance_fluid_24hr} mL + {burn_fluid_24hr} mL = {total_fluid_24hr} mL")
    print("\nFinal Calculation:")
    print(f"Average Hourly Rate = {total_fluid_24hr} mL / 24 hours")
    print(f"Result = {average_hourly_rate:.2f} cc/hr")

calculate_pediatric_burn_fluid()