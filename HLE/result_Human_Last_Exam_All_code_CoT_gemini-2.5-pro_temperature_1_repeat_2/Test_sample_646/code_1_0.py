import sys

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming contaminated produce.
    The function follows these steps:
    1. Calculates the initial concentration of PFHxS in the soil.
    2. Calculates the concentration of PFHxS in the soil solution.
    3. Calculates the concentration of PFHxS in the consumed fruits and legumes.
    4. Calculates the total Chronic Daily Intake (CDI) from both food sources.
    5. Calculates the final Hazard Quotient (HQ) by comparing the CDI to the reference dose.
    """

    # --- Given Data ---
    # Soil and Contamination
    area_m2 = 250000
    depth_m = 0.6
    foc = 0.03  # fraction of organic carbon
    theta_w = 0.35  # volumetric water content (L/L)
    bulk_density_kg_m3 = 1500
    foam_volume_L = 1000
    foam_conc_ug_L = 1000000
    koc_L_kg = 500 # Assumed typical organic carbon-water partition coefficient for PFHxS

    # Exposure Data (Male)
    body_weight_kg = 80
    rfd_ug_kg_day = 0.02  # Reference Dose

    # Fruit Consumption
    ir_fruit_g_day = 300
    bf_fruit = 0.5
    puf_fruit = 0.1
    tscf_fruit = 5

    # Legume Consumption
    ir_legume_g_day = 50
    bf_legume = 0.3
    puf_legume = 0.2
    tscf_legume = 5

    # --- Step 1 & 2: Calculate PFHxS Concentration in Soil (C_soil) ---
    print("Step 1 & 2: Calculating PFHxS Concentration in Soil (C_soil)")
    
    total_pfhxs_ug = foam_volume_L * foam_conc_ug_L
    soil_volume_m3 = area_m2 * depth_m
    soil_mass_kg = soil_volume_m3 * bulk_density_kg_m3
    c_soil_ug_kg = total_pfhxs_ug / soil_mass_kg
    
    print(f"Total PFHxS applied = {foam_volume_L:.0f} L * {foam_conc_ug_L:,.0f} ug/L = {total_pfhxs_ug:,.0f} ug")
    print(f"Total soil mass = {soil_volume_m3:,.0f} m^3 * {bulk_density_kg_m3:.0f} kg/m^3 = {soil_mass_kg:,.0f} kg")
    print(f"Concentration in soil (C_soil) = {total_pfhxs_ug:,.0f} ug / {soil_mass_kg:,.0f} kg = {c_soil_ug_kg:.4f} ug/kg\n")

    # --- Step 3: Calculate PFHxS Concentration in Soil Solution (C_solution) ---
    print("Step 3: Calculating PFHxS Concentration in Soil Solution (C_solution)")
    
    kd_L_kg = koc_L_kg * foc
    bulk_density_kg_L = bulk_density_kg_m3 / 1000.0
    c_solution_ug_L = c_soil_ug_kg / (kd_L_kg + (theta_w / bulk_density_kg_L))
    
    print(f"Soil-water partition coefficient (Kd) = {koc_L_kg:.0f} L/kg (Koc) * {foc} (foc) = {kd_L_kg:.2f} L/kg")
    print(f"Concentration in soil solution (C_solution) = {c_soil_ug_kg:.4f} ug/kg / ({kd_L_kg:.2f} L/kg + ({theta_w} / {bulk_density_kg_L} kg/L)) = {c_solution_ug_L:.4f} ug/L\n")

    # --- Step 4: Calculate PFHxS Concentration in Produce (C_plant) ---
    print("Step 4: Calculating PFHxS Concentration in Produce (assuming 1L plant tissue ≈ 1kg)")
    
    c_fruit_ug_kg = c_solution_ug_L * tscf_fruit * puf_fruit
    c_legume_ug_kg = c_solution_ug_L * tscf_legume * puf_legume
    
    print(f"Concentration in fruit = {c_solution_ug_L:.4f} ug/L * {tscf_fruit} (TSCF) * {puf_fruit} (PUF) = {c_fruit_ug_kg:.4f} ug/kg")
    print(f"Concentration in legume = {c_solution_ug_L:.4f} ug/L * {tscf_legume} (TSCF) * {puf_legume} (PUF) = {c_legume_ug_kg:.4f} ug/kg\n")
    
    # --- Step 5: Calculate Chronic Daily Intake (CDI) ---
    print("Step 5: Calculating Chronic Daily Intake (CDI)")
    
    ir_fruit_kg_day = ir_fruit_g_day / 1000.0
    cdi_fruit = (c_fruit_ug_kg * ir_fruit_kg_day * bf_fruit) / body_weight_kg
    
    ir_legume_kg_day = ir_legume_g_day / 1000.0
    cdi_legume = (c_legume_ug_kg * ir_legume_kg_day * bf_legume) / body_weight_kg
    
    cdi_total = cdi_fruit + cdi_legume
    
    print(f"CDI from fruit = ({c_fruit_ug_kg:.4f} ug/kg * {ir_fruit_kg_day} kg/day * {bf_fruit}) / {body_weight_kg} kg = {cdi_fruit:.8f} ug/kg/day")
    print(f"CDI from legume = ({c_legume_ug_kg:.4f} ug/kg * {ir_legume_kg_day} kg/day * {bf_legume}) / {body_weight_kg} kg = {cdi_legume:.8f} ug/kg/day")
    print(f"Total CDI = {cdi_fruit:.8f} + {cdi_legume:.8f} = {cdi_total:.8f} ug/kg/day\n")

    # --- Step 6: Calculate Hazard Quotient (HQ) ---
    print("Step 6: Calculating the Final Hazard Quotient (HQ)")
    
    hq = cdi_total / rfd_ug_kg_day
    
    print(f"Hazard Quotient (HQ) = Total CDI / RfD")
    print(f"HQ = {cdi_total:.8f} ug/kg/day / {rfd_ug_kg_day} ug/kg/day = {hq:.4f}\n")
    
    # Store the final answer in a specific format for parsing
    sys.stdout.write(f'<<<>>>\n') # Guard for answer extraction
    sys.stdout.write(f'<<<{hq}>>>\n')

calculate_hazard_quotient()
