import sys

def calculate_hazard_quotient():
    """
    This function calculates the hazard quotient for a man exposed to PFHxS
    through consumption of contaminated produce.
    """
    # --- Step 1: Define all constants from the problem ---

    # Contamination source
    volume_foam_L = 1000.0
    conc_in_foam_ug_per_L = 1000000.0

    # Soil characteristics
    area_m2 = 250000.0
    depth_m = 0.6
    bulk_density_kg_per_m3 = 1500.0

    # Human exposure characteristics
    body_weight_kg = 80.0
    ref_dose_ug_per_kg_day = 0.02

    # Fruit consumption details
    intake_fruit_g_day = 300.0
    baf_fruit = 0.5
    puf_fruit = 0.1

    # Legume consumption details
    intake_legume_g_day = 50.0
    baf_legume = 0.3
    puf_legume = 0.2

    # --- Step 2: Intermediate Calculations ---

    # Calculate total mass of PFHxS applied to the soil
    total_pfhxs_ug = volume_foam_L * conc_in_foam_ug_per_L

    # Calculate the total mass of the contaminated soil
    volume_soil_m3 = area_m2 * depth_m
    mass_soil_kg = volume_soil_m3 * bulk_density_kg_per_m3

    # Calculate the concentration of PFHxS in the soil (C_soil)
    c_soil_ug_per_kg = total_pfhxs_ug / mass_soil_kg

    # Calculate the concentration of PFHxS in the produce
    c_fruit_ug_per_kg = c_soil_ug_per_kg * puf_fruit
    c_legume_ug_per_kg = c_soil_ug_per_kg * puf_legume

    # Convert daily intake from grams to kilograms
    intake_fruit_kg_day = intake_fruit_g_day / 1000.0
    intake_legume_kg_day = intake_legume_g_day / 1000.0

    # Calculate the total absorbed dose of PFHxS per day
    absorbed_dose_fruit_ug_day = c_fruit_ug_per_kg * intake_fruit_kg_day * baf_fruit
    absorbed_dose_legume_ug_day = c_legume_ug_per_kg * intake_legume_kg_day * baf_legume
    total_absorbed_dose_ug_day = absorbed_dose_fruit_ug_day + absorbed_dose_legume_ug_day

    # Calculate the Systemic Exposure Dose (SED)
    sed_ug_per_kg_day = total_absorbed_dose_ug_day / body_weight_kg

    # Calculate the final Hazard Quotient (HQ)
    hazard_quotient = sed_ug_per_kg_day / ref_dose_ug_per_kg_day

    # --- Step 3: Print the calculation steps and the final answer ---

    print("--- Calculation of PFHxS Concentration in Soil (C_soil) ---")
    print(f"Total PFHxS applied = {volume_foam_L} L * {conc_in_foam_ug_per_L} µg/L = {total_pfhxs_ug:,.0f} µg")
    print(f"Total soil mass = {area_m2:,.0f} m² * {depth_m} m * {bulk_density_kg_per_m3} kg/m³ = {mass_soil_kg:,.0f} kg")
    print(f"C_soil = {total_pfhxs_ug:,.0f} µg / {mass_soil_kg:,.0f} kg = {c_soil_ug_per_kg:.4f} µg/kg\n")
    
    print("--- Calculation of Systemic Exposure Dose (SED) ---")
    print("Daily absorbed dose = (C_fruit * IR_fruit * BAF_fruit) + (C_legume * IR_legume * BAF_legume)")
    print(f"Daily absorbed dose = (({c_soil_ug_per_kg:.4f} * {puf_fruit}) * {intake_fruit_kg_day} * {baf_fruit}) + (({c_soil_ug_per_kg:.4f} * {puf_legume}) * {intake_legume_kg_day} * {baf_legume})")
    print(f"Daily absorbed dose = {absorbed_dose_fruit_ug_day:.4f} µg/day (from fruits) + {absorbed_dose_legume_ug_day:.4f} µg/day (from legumes) = {total_absorbed_dose_ug_day:.4f} µg/day")
    print(f"SED = Total daily absorbed dose / Body weight")
    print(f"SED = {total_absorbed_dose_ug_day:.4f} µg/day / {body_weight_kg} kg = {sed_ug_per_kg_day:.6f} µg/kg/day\n")

    print("--- Final Hazard Quotient (HQ) Calculation ---")
    print(f"HQ = Systemic Exposure Dose (SED) / Reference Dose (RfD)")
    print(f"HQ = {sed_ug_per_kg_day:.6f} µg/kg/day / {ref_dose_ug_per_kg_day} µg/kg/day")
    print(f"The final Hazard Quotient is: {hazard_quotient:.2f}")

    # The final answer in the required format
    sys.stdout.write(f"\n<<<{hazard_quotient:.2f}>>>")

if __name__ == '__main__':
    calculate_hazard_quotient()
