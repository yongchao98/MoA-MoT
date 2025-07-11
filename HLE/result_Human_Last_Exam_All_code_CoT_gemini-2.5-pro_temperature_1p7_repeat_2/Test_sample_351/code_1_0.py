import math

def calculate_pediatric_burn_fluids():
    """
    Calculates the maintenance fluid needs for a pediatric burn patient
    using the Galveston formula.
    """
    # Patient data
    weight_kg = 25
    tbsa_percent = 45
    bsa_m2 = 0.95

    # --- Calculation Plan ---
    print("This script calculates the total hourly fluid rate for a pediatric burn patient using the Galveston formula.")
    print("The Galveston formula provides a single, constant infusion rate for the first 24 hours post-burn.")
    print(f"Patient Data: 7 years old, {weight_kg}kg, {tbsa_percent}% TBSA burn, {bsa_m2} m^2 BSA.\n")

    # Step 1: Calculate Burned Body Surface Area (BSA_burn)
    print("Step 1: Calculate the Burned Body Surface Area (BSA_burn) in m^2.")
    tbsa_decimal = tbsa_percent / 100
    bsa_burn = bsa_m2 * tbsa_decimal
    print(f"   BSA_burn = Total BSA * (%TBSA / 100) = {bsa_m2} m^2 * {tbsa_decimal} = {bsa_burn:.4f} m^2\n")

    # Step 2: Calculate the 24-hour fluid volume for the burn component
    print("Step 2: Calculate the 24-hour fluid for the burn area (5000 mL/m^2).")
    fluid_for_burn = 5000 * bsa_burn
    print(f"   Fluid for Burn = 5000 mL/m^2 * BSA_burn = 5000 * {bsa_burn:.4f} = {fluid_for_burn:.2f} mL\n")

    # Step 3: Calculate the 24-hour fluid volume for the maintenance component
    print("Step 3: Calculate the 24-hour maintenance fluid (2000 mL/m^2).")
    fluid_for_maintenance = 2000 * bsa_m2
    print(f"   Fluid for Maintenance = 2000 mL/m^2 * Total BSA = 2000 * {bsa_m2} = {fluid_for_maintenance:.2f} mL\n")

    # Step 4: Calculate the total 24-hour fluid volume
    print("Step 4: Sum the components to find the total 24-hour fluid requirement.")
    total_fluid_24hr = fluid_for_burn + fluid_for_maintenance
    print(f"   Total 24hr Fluid = {fluid_for_burn:.2f} mL + {fluid_for_maintenance:.2f} mL = {total_fluid_24hr:.2f} mL\n")

    # Step 5: Calculate the final hourly rate
    print("Step 5: Divide the total 24-hour volume by 24 to get the hourly rate in cc/hr.")
    hourly_rate = total_fluid_24hr / 24
    print(f"   Hourly Rate = {total_fluid_24hr:.2f} mL / 24 hr = {hourly_rate:.2f} cc/hr\n")
    
    # Final equation summary with all numbers
    print("--- Final Equation Summary ---")
    print(f"Hourly Rate = ((5000 * ({bsa_m2} * {tbsa_decimal})) + (2000 * {bsa_m2})) / 24")
    print(f"Result: {hourly_rate:.2f} cc/hr")

calculate_pediatric_burn_fluids()