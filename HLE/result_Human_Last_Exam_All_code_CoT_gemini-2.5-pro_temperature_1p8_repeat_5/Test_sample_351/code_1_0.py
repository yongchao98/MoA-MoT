import math

def calculate_pediatric_burn_fluids():
    """
    Calculates and displays the maintenance fluid needs for a pediatric burn patient.

    This calculation uses the Body Surface Area (BSA) method, which accounts for
    both basal fluid needs and evaporative losses from the burn.
    """
    # Patient data from the prompt
    weight_kg = 25
    tbsa_percent = 45
    bsa_m2 = 0.95

    # Constants for the BSA method formula
    BASAL_RATE_PER_M2 = 1500  # mL/m^2/day
    EVAPORATIVE_RATE_PER_M2_BURN = 3500  # mL/m^2 of burn area/day
    HOURS_IN_DAY = 24

    # --- Calculation Steps ---

    # 1. Calculate basal fluid needs based on total BSA
    basal_fluid_per_day = BASAL_RATE_PER_M2 * bsa_m2

    # 2. Calculate the size of the burn area in square meters
    tbsa_decimal = tbsa_percent / 100
    burn_area_m2 = bsa_m2 * tbsa_decimal

    # 3. Calculate additional fluid needed for evaporative loss from the burn area
    evaporative_fluid_per_day = EVAPORATIVE_RATE_PER_M2_BURN * burn_area_m2

    # 4. Calculate the total fluid needed over 24 hours
    total_fluid_per_day = basal_fluid_per_day + evaporative_fluid_per_day

    # 5. Convert the 24-hour total into an hourly rate
    hourly_rate = total_fluid_per_day / HOURS_IN_DAY

    # --- Output Results ---
    print("Calculating Total Fluid Needs for a Pediatric Burn Patient (BSA Method)")
    print(f"Patient Data: {weight_kg}kg, {tbsa_percent}% TBSA burn, {bsa_m2} m^2 BSA")
    print("-" * 60)
    print("Step 1: Calculate Basal Fluid Needs (mL/day)")
    print(f"   {BASAL_RATE_PER_M2} mL/m^2 * {bsa_m2} m^2 = {basal_fluid_per_day:.2f} mL/day")
    print("-" * 60)
    print("Step 2: Calculate Evaporative Fluid Loss (mL/day)")
    print(f"   Burn Area = {bsa_m2} m^2 (BSA) * {tbsa_decimal} ({tbsa_percent}%) = {burn_area_m2:.4f} m^2")
    print(f"   {EVAPORATIVE_RATE_PER_M2_BURN} mL/m^2 * {burn_area_m2:.4f} m^2 (Burn Area) = {evaporative_fluid_per_day:.2f} mL/day")
    print("-" * 60)
    print("Step 3: Calculate Total Daily Fluid Needs (mL/day)")
    print(f"   {basal_fluid_per_day:.2f} (Basal) + {evaporative_fluid_per_day:.2f} (Evaporative) = {total_fluid_per_day:.2f} mL/day")
    print("-" * 60)
    print("Step 4: Calculate Final Hourly Rate (cc/hr)")
    print("   The complete formula is: ((Basal Rate * BSA) + (Evaporative Rate * Burn Area)) / 24")
    print(f"   (({BASAL_RATE_PER_M2} * {bsa_m2}) + ({EVAPORATIVE_RATE_PER_M2_BURN} * ({bsa_m2} * {tbsa_decimal}))) / {HOURS_IN_DAY}")
    print(f"\n   The final required fluid rate is {hourly_rate:.2f} cc/hr.")


calculate_pediatric_burn_fluids()
<<<121.72>>>