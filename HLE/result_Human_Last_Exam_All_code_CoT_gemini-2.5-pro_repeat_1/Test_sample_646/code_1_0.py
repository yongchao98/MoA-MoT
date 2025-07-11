import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming contaminated produce.
    """
    # Step 1: Define given parameters

    # Site and Soil Characteristics
    area_m2 = 250000.0  # m²
    depth_m = 0.6  # m
    f_oc = 0.03  # organic carbon content (3%)
    theta_w = 0.35  # volumetric water content (L water / L soil)
    rho_b_kg_m3 = 1500.0  # soil bulk density in kg/m³

    # Contamination Details
    volume_foam_L = 1000.0  # L
    conc_foam_ug_L = 1000000.0  # μg/L

    # Exposure Scenario
    body_weight_kg = 80.0  # kg
    rfd = 0.02  # Reference Dose in µg/kg/day

    # Fruit Consumption
    intake_fruit_g_day = 300.0  # g/day
    bf_fruit = 0.5  # bioavailability factor
    puf_fruit = 0.1  # plant uptake factor
    tscf_fruit = 5.0  # transpiration stream concentration factor

    # Legume Consumption
    intake_legume_g_day = 50.0  # g/day
    bf_legume = 0.3  # bioavailability factor
    puf_legume = 0.2  # plant uptake factor
    tscf_legume = 5.0  # transpiration stream concentration factor

    # Assumption for PFHxS chemical property
    # The octanol-carbon partition coefficient (Koc) is not given.
    # A literature value for log Koc of PFHxS is approximately 2.8.
    log_koc_pfhxs = 2.8
    koc = 10**log_koc_pfhxs  # L/kg

    # --- Calculations ---

    # Calculate total mass of PFHxS applied
    total_mass_pfhxs_ug = conc_foam_ug_L * volume_foam_L

    # Calculate volume and mass of contaminated soil
    volume_soil_m3 = area_m2 * depth_m
    mass_soil_kg = volume_soil_m3 * rho_b_kg_m3

    # Calculate concentration of PFHxS in soil (C_soil)
    c_soil_ug_kg = total_mass_pfhxs_ug / mass_soil_kg

    # Calculate soil-water partition coefficient (Kd)
    kd_L_kg = koc * f_oc

    # Calculate concentration of PFHxS in soil solution (C_solution)
    # Convert bulk density to kg/L: 1 m³ = 1000 L
    rho_b_kg_L = rho_b_kg_m3 / 1000.0
    # Using the formula: C_solution (ug/L) = C_soil (ug/kg) / (Kd (L/kg) + theta_w / rho_b (kg/L))
    c_solution_ug_L = c_soil_ug_kg / (kd_L_kg + (theta_w / rho_b_kg_L))

    # Calculate PFHxS concentration in produce
    # Assuming 1 L of water uptake corresponds to 1 kg of fresh plant weight
    c_fruit_ug_kg = c_solution_ug_L * tscf_fruit * puf_fruit
    c_legume_ug_kg = c_solution_ug_L * tscf_legume * puf_legume

    # Calculate average daily dose (ADD) from each food source
    # Convert daily intake from g/day to kg/day
    intake_fruit_kg_day = intake_fruit_g_day / 1000.0
    intake_legume_kg_day = intake_legume_g_day / 1000.0
    # Formula: ADD (ug/kg/day) = (C_plant * Intake_food * BF) / BW
    add_fruit = (c_fruit_ug_kg * intake_fruit_kg_day * bf_fruit) / body_weight_kg
    add_legume = (c_legume_ug_kg * intake_legume_kg_day * bf_legume) / body_weight_kg

    # Calculate total average daily dose (ADD_total)
    add_total = add_fruit + add_legume

    # Calculate Hazard Quotient (HQ)
    # Formula: HQ = ADD_total / RfD
    hq = add_total / rfd

    # --- Output Results ---
    print("This calculation determines the Hazard Quotient (HQ) for an individual exposed to PFHxS through consumption of contaminated produce.")
    print("An assumption for the PFHxS soil-organic carbon partition coefficient (log Koc = 2.8) was made based on literature values.\n")
    print("--- Final Hazard Quotient Calculation ---\n")
    
    # Print the breakdown of the total daily intake
    print(f"Daily Dose from Fruits = (Concentration in Fruit * Daily Intake * Bioavailability) / Body Weight")
    print(f"Daily Dose from Fruits = ({c_fruit_ug_kg:.4f} µg/kg * {intake_fruit_kg_day} kg/day * {bf_fruit}) / {body_weight_kg} kg = {add_fruit:.8f} µg/kg/day\n")
    
    print(f"Daily Dose from Legumes = (Concentration in Legumes * Daily Intake * Bioavailability) / Body Weight")
    print(f"Daily Dose from Legumes = ({c_legume_ug_kg:.4f} µg/kg * {intake_legume_kg_day} kg/day * {bf_legume}) / {body_weight_kg} kg = {add_legume:.8f} µg/kg/day\n")
    
    # Print the final equation with all numbers
    print(f"Total Average Daily Dose (ADD) = Dose from Fruits + Dose from Legumes")
    print(f"Total Average Daily Dose (ADD) = {add_fruit:.8f} + {add_legume:.8f} = {add_total:.8f} µg/kg/day\n")
    
    print(f"Hazard Quotient (HQ) = Total Average Daily Dose / Reference Dose\n")
    print(f"HQ = {add_total:.8f} / {rfd} = {hq:.4f}")
    
    return hq

# Execute the calculation and store the final answer
final_hq = calculate_hazard_quotient()
print(f"\n<<<{final_hq:.4f}>>>")