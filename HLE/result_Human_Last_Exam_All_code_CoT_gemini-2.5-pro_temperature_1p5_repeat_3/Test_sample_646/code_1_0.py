import math

def calculate_hazard_quotient():
    """
    Calculates the hazard quotient for a man exposed to PFHxS from contaminated produce.
    """
    # --- Input Data ---
    # Soil and Site characteristics
    area_m2 = 250000
    depth_m = 0.6
    f_oc = 0.03  # 3% organic carbon
    theta_w = 0.35  # L water/L soil
    bulk_density_kg_m3 = 1500

    # Contamination details
    foam_volume_L = 1000
    pfhxs_conc_in_foam_ug_L = 1000000
    log_koc_pfhxs = 3.8 # Standard literature value for log Koc of PFHxS

    # Exposure scenario
    bw_kg = 80 # body weight in kg

    # Fruits consumption
    ir_fruits_g_day = 300
    bf_fruits = 0.5
    puf_fruits = 0.1
    tscf_fruits = 5

    # Legumes consumption
    ir_legumes_g_day = 50
    bf_legumes = 0.3
    puf_legumes = 0.2
    tscf_legumes = 5

    # Toxicological reference value
    rfd_ug_kg_day = 0.02

    # --- Step 1: Calculate total mass of contaminant ---
    total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_L

    # --- Step 2: Calculate total mass of soil ---
    soil_volume_m3 = area_m2 * depth_m
    soil_mass_kg = soil_volume_m3 * bulk_density_kg_m3

    # --- Step 3: Calculate PFHxS concentration in soil (C_soil) ---
    c_soil_ug_kg = total_pfhxs_mass_ug / soil_mass_kg

    # --- Step 4: Calculate PFHxS concentration in soil water (C_w) ---
    # Calculate K_oc from log K_oc
    k_oc_L_kg = 10**log_koc_pfhxs
    # Calculate soil-water partition coefficient (K_d)
    k_d_L_kg = k_oc_L_kg * f_oc
    # Convert bulk density to kg/L for unit consistency in the formula
    bulk_density_kg_L = bulk_density_kg_m3 / 1000.0
    # C_w = C_soil / (K_d + θ_w / ρ_b)
    c_w_ug_L = c_soil_ug_kg / (k_d_L_kg + (theta_w / bulk_density_kg_L))

    # --- Step 5: Calculate PFHxS concentration in produce ---
    # Concentration in fruits (assuming fresh weight kg is approx. equal to L)
    c_fruits_ug_kg = c_w_ug_L * tscf_fruits * puf_fruits
    # Concentration in legumes (assuming fresh weight kg is approx. equal to L)
    c_legumes_ug_kg = c_w_ug_L * tscf_legumes * puf_legumes

    # --- Step 6: Calculate Average Daily Dose (ADD) ---
    # Convert intake rates from g/day to kg/day
    ir_fruits_kg_day = ir_fruits_g_day / 1000.0
    ir_legumes_kg_day = ir_legumes_g_day / 1000.0
    
    # Absorbed dose from each food type
    dose_fruits_ug_day = c_fruits_ug_kg * ir_fruits_kg_day * bf_fruits
    dose_legumes_ug_day = c_legumes_ug_kg * ir_legumes_kg_day * bf_legumes
    
    # Total absorbed dose divided by body weight
    add_ug_kg_day = (dose_fruits_ug_day + dose_legumes_ug_day) / bw_kg

    # --- Step 7: Calculate Hazard Quotient (HQ) ---
    hq = add_ug_kg_day / rfd_ug_kg_day

    # --- Print Results ---
    print("--- Hazard Quotient Calculation ---")
    print("\n1. Average Daily Dose (ADD) Calculation:")
    print(f"The Average Daily Dose is the sum of absorbed doses from all food sources, divided by body weight.")
    print(f"ADD = [ (C_fruits * IR_fruits * BF_fruits) + (C_legumes * IR_legumes * BF_legumes) ] / BW")
    print(f"ADD = [ ({c_fruits_ug_kg:.5f} ug/kg * {ir_fruits_kg_day} kg/day * {bf_fruits}) + ({c_legumes_ug_kg:.5f} ug/kg * {ir_legumes_kg_day} kg/day * {bf_legumes}) ] / {bw_kg} kg")
    print(f"ADD = {add_ug_kg_day:.8f} µg/kg/day")

    print("\n2. Hazard Quotient (HQ) Calculation:")
    print("The Hazard Quotient is the Average Daily Dose divided by the Reference Dose.")
    print(f"HQ = ADD / RfD")
    print(f"HQ = {add_ug_kg_day:.8f} / {rfd_ug_kg_day}")
    
    # Final Answer
    print(f"\nThe Hazard Quotient is {hq:.4f}")
    print(f"<<<{hq:.4f}>>>")

calculate_hazard_quotient()