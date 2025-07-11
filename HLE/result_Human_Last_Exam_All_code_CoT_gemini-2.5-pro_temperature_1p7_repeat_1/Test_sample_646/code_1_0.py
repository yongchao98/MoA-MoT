import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming contaminated produce.
    """
    # --- 1. Given Data ---
    # Contamination Data
    volume_foam_L = 1000  # L
    conc_foam_ug_per_L = 1000000  # μg/L
    area_m2 = 250000  # m²
    depth_m = 0.6  # m

    # Soil Properties
    f_oc = 0.03  # organic carbon content (3%)
    theta_w = 0.35  # volumetric water content (L water/L soil)
    rho_b_kg_per_m3 = 1500  # bulk density (kg/m³)

    # Exposure Scenario
    bw_kg = 80  # kg

    # Consumption Data (Fruits)
    ir_fruit_g_per_day = 300  # g/day
    bf_fruit = 0.5
    tscf_fruit = 5.0

    # Consumption Data (Legumes)
    ir_legume_g_per_day = 50  # g/day
    bf_legume = 0.3
    tscf_legume = 5.0

    # Toxicological Data
    rfd_ug_per_kg_day = 0.02  # μg/kg/day
    
    # Assumed value: log K_oc for PFHxS, as it's not provided in the problem.
    # This is a common practice in risk assessment when data is missing.
    log_koc_pfhxs = 3.2
    
    # --- 2. Calculations ---
    # Total mass of PFHxS applied
    mass_pfhxs_ug = volume_foam_L * conc_foam_ug_per_L
    
    # Total volume and mass of contaminated soil
    volume_soil_m3 = area_m2 * depth_m
    mass_soil_kg = volume_soil_m3 * rho_b_kg_per_m3
    
    # Concentration in soil (Cs)
    c_s_ug_per_kg = mass_pfhxs_ug / mass_soil_kg
    
    # Concentration in soil water (Cw)
    k_oc = 10**log_koc_pfhxs  # L/kg
    k_d = k_oc * f_oc  # L/kg
    # Convert bulk density to kg/L for unit consistency
    rho_b_kg_per_L = rho_b_kg_per_m3 / 1000
    # Cw = Cs / (Kd + θw/ρb)
    c_w_ug_per_L = c_s_ug_per_kg / (k_d + (theta_w / rho_b_kg_per_L))
    
    # Concentration in produce (C_plant)
    # Assuming C_plant_fw (μg/kg) ≈ Cw (μg/L) * TSCF, as 1L of plant matter ≈ 1kg
    c_fruit_ug_per_kg = c_w_ug_per_L * tscf_fruit
    c_legume_ug_per_kg = c_w_ug_per_L * tscf_legume
    
    # Daily Intake from each food source
    # Convert intake rates from g/day to kg/day
    ir_fruit_kg_per_day = ir_fruit_g_per_day / 1000
    ir_legume_kg_per_day = ir_legume_g_per_day / 1000
    
    intake_fruit_ug_day = c_fruit_ug_per_kg * ir_fruit_kg_per_day * bf_fruit
    intake_legume_ug_day = c_legume_ug_per_kg * ir_legume_kg_per_day * bf_legume
    
    # Total daily intake (absorbed dose)
    total_daily_intake_ug_day = intake_fruit_ug_day + intake_legume_ug_day
    
    # Chronic Daily Intake (CDI) per kg of body weight
    cdi_ug_per_kg_day = total_daily_intake_ug_day / bw_kg
    
    # Hazard Quotient (HQ)
    hq = cdi_ug_per_kg_day / rfd_ug_per_kg_day

    # --- 3. Print Results ---
    print("--- Calculation Steps ---")
    print(f"1. Total mass of PFHxS contaminant: {mass_pfhxs_ug:,.0f} µg")
    print(f"2. Total mass of contaminated soil: {mass_soil_kg:,.0f} kg")
    print(f"3. Concentration in soil (Cs): {c_s_ug_per_kg:.4f} µg/kg")
    print(f"4. Concentration in soil water (Cw): {c_w_ug_per_L:.4f} µg/L")
    print(f"5. Concentration in fruit (C_fruit): {c_fruit_ug_per_kg:.4f} µg/kg")
    print(f"   Concentration in legume (C_legume): {c_legume_ug_per_kg:.4f} µg/kg")
    print(f"6. Total daily intake (absorbed): {total_daily_intake_ug_day:.4f} µg/day")
    print(f"7. Chronic Daily Intake (CDI): {cdi_ug_per_kg_day:.6f} µg/kg/day")
    print("\n--- Final Hazard Quotient Calculation ---")
    print("HQ = CDI / RfD")
    print("HQ = (Total Daily Intake / Body Weight) / Reference Dose")
    final_eq_str = (
        f"HQ = (({c_fruit_ug_per_kg:.4f} * {ir_fruit_kg_per_day} * {bf_fruit}) + "
        f"({c_legume_ug_per_kg:.4f} * {ir_legume_kg_per_day} * {bf_legume})) / "
        f"{bw_kg} / {rfd_ug_per_kg_day}"
    )
    print(final_eq_str)
    print(f"HQ = {hq:.4f}")
    
    return hq

# Execute the function
calculated_hq = calculate_hazard_quotient()
# Final answer in the required format
print(f"<<<{calculated_hq:.4f}>>>")
