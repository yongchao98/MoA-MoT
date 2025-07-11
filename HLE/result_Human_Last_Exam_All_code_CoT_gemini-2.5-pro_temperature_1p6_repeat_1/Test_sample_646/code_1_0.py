import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming contaminated produce.
    """
    # Step 1: Define constants and inputs from the problem description
    # Soil and contaminant properties
    area_m2 = 250000
    depth_m = 0.6
    organic_carbon_content_percent = 3.0
    volumetric_water_content_ratio = 0.35  # L water / L soil
    bulk_density_kg_m3 = 1500
    foam_volume_L = 1000
    pfhxs_in_foam_ug_L = 1000000

    # Exposure scenario properties
    body_weight_kg = 80.0
    ref_dose_ug_kg_day = 0.02

    # Fruit properties
    fruit_intake_g_day = 300.0
    fruit_bioavailability = 0.5
    fruit_puf = 0.1
    fruit_tscf = 5.0

    # Legume properties
    legume_intake_g_day = 50.0
    legume_bioavailability = 0.3
    legume_puf = 0.2
    legume_tscf = 5.0

    # ASSUMPTION: The soil organic carbon-water partitioning coefficient (Koc) for PFHxS
    # is not provided. A literature-based log Koc value of 3.5 is assumed.
    log_koc_pfhxs = 3.5
    koc_L_kg = 10**log_koc_pfhxs

    # Step 2: Calculate initial concentration in soil (Cs)
    total_pfhxs_ug = foam_volume_L * pfhxs_in_foam_ug_L
    soil_volume_m3 = area_m2 * depth_m
    soil_mass_kg = soil_volume_m3 * bulk_density_kg_m3
    cs_ug_kg = total_pfhxs_ug / soil_mass_kg

    # Step 3: Calculate concentration in soil solution (Cw)
    foc = organic_carbon_content_percent / 100.0
    kd_L_kg = koc_L_kg * foc
    # Convert bulk density to kg/L
    bulk_density_kg_L = bulk_density_kg_m3 / 1000.0
    # Convert volumetric water content (theta_w) to L/kg to match Kd units
    theta_w_L_kg = volumetric_water_content_ratio / bulk_density_kg_L
    # Concentration in soil water (Cw) in ug/L
    cw_ug_L = cs_ug_kg / (kd_L_kg + theta_w_L_kg)

    # Step 4: Calculate plant concentrations (µg/g)
    # Concentration in fruit, assuming 1L plant tissue ≈ 1kg, then converting kg to g
    c_fruit_ug_g = (cw_ug_L * fruit_tscf * fruit_puf) / 1000.0
    # Concentration in legumes
    c_legume_ug_g = (cw_ug_L * legume_tscf * legume_puf) / 1000.0

    # Step 5: Calculate total absorbed daily intake (DI)
    # Absorbed daily dose from fruit and legumes in µg/day
    dose_fruit_ug_day = c_fruit_ug_g * fruit_intake_g_day * fruit_bioavailability
    dose_legume_ug_day = c_legume_ug_g * legume_intake_g_day * legume_bioavailability
    
    total_absorbed_dose_ug_day = dose_fruit_ug_day + dose_legume_ug_day

    # Daily Intake (DI) per kg body weight
    daily_intake_ug_kg_day = total_absorbed_dose_ug_day / body_weight_kg

    # Step 6: Calculate the Hazard Quotient (HQ)
    hazard_quotient = daily_intake_ug_kg_day / ref_dose_ug_kg_day

    # Step 7: Print the final equations and result
    print("Calculating Daily Intake (DI):")
    print(f"DI = (Absorbed Dose from Fruit + Absorbed Dose from Legumes) / Body Weight")
    print(f"Absorbed Dose from Fruit = Concentration [µg/g] * Intake [g/day] * Bioavailability")
    print(f"                       = {c_fruit_ug_g:.8f} * {fruit_intake_g_day} * {fruit_bioavailability} = {dose_fruit_ug_day:.8f} µg/day")
    print(f"Absorbed Dose from Legumes = {c_legume_ug_g:.8f} * {legume_intake_g_day} * {legume_bioavailability} = {dose_legume_ug_day:.8f} µg/day")
    print(f"Total Absorbed Dose = {dose_fruit_ug_day:.8f} + {dose_legume_ug_day:.8f} = {total_absorbed_dose_ug_day:.8f} µg/day")
    print(f"DI = {total_absorbed_dose_ug_day:.8f} µg/day / {body_weight_kg} kg = {daily_intake_ug_kg_day:.8f} µg/kg/day\n")

    print("Calculating Hazard Quotient (HQ):")
    print("HQ = Daily Intake (DI) / Reference Dose (RfD)")
    print(f"HQ = {daily_intake_ug_kg_day:.8f} µg/kg/day / {ref_dose_ug_kg_day} µg/kg/day")
    print(f"Final Hazard Quotient = {hazard_quotient:.4f}")
    
    # Returning the final numerical answer for the platform.
    return hazard_quotient

# Execute the calculation and print the result.
calculated_hq = calculate_hazard_quotient()
# The final answer is formatted as requested below.
# <<<final answer>>>
final_answer = f"<<<{calculated_hq:.4f}>>>"
# This print is for display in a standard terminal.
# The 'final_answer' variable holds the required format.
# print(final_answer)

if __name__ == '__main__':
    # This block is for standard execution and testing.
    # The platform extracts the answer from the <<<>>> format.
    pass
