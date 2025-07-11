import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient for PFHxS exposure from consuming produce
    grown on contaminated land.
    """
    # --- Step 0: Define constants and given data ---
    # Soil and Site Data
    area_m2 = 250000
    depth_m = 0.6
    f_oc = 0.03  # organic carbon content fraction
    theta_w = 0.35  # volumetric water content (L water/L soil)
    rho_b_kg_m3 = 1500  # bulk density in kg/m^3

    # Contaminant Data
    foam_volume_L = 1000
    foam_conc_ug_L = 1000000
    log_koc_PFHxS = 3.2 # Standard literature value for log Koc of PFHxS
    koc = 10**log_koc_PFHxS

    # Exposure Data
    body_weight_kg = 80
    rfd_ug_kg_day = 0.02 # Reference Dose

    # Fruits Data
    ir_fruit_g_day = 300
    bf_fruit = 0.5
    puf_fruit = 0.1
    tscf = 5 # Transpiration Stream Concentration Factor (same for both)

    # Legumes Data
    ir_legume_g_day = 50
    bf_legume = 0.3
    puf_legume = 0.2

    # --- Unit Conversions ---
    ir_fruit_kg_day = ir_fruit_g_day / 1000
    ir_legume_kg_day = ir_legume_g_day / 1000
    rho_b_kg_L = rho_b_kg_m3 / 1000 # Convert bulk density to kg/L

    print("--- Step 1: Calculate PFHxS Concentration in Soil (Cs) ---")
    total_mass_pfhxs_ug = foam_volume_L * foam_conc_ug_L
    print(f"Total mass of PFHxS applied = {foam_volume_L} L * {foam_conc_ug_L} µg/L = {total_mass_pfhxs_ug:e} µg")

    soil_volume_m3 = area_m2 * depth_m
    soil_mass_kg = soil_volume_m3 * rho_b_kg_m3
    print(f"Total mass of affected soil = ({area_m2} m² * {depth_m} m) * {rho_b_kg_m3} kg/m³ = {soil_mass_kg:e} kg")

    cs_ug_kg = total_mass_pfhxs_ug / soil_mass_kg
    print(f"Concentration in soil (Cs) = {total_mass_pfhxs_ug:e} µg / {soil_mass_kg:e} kg = {cs_ug_kg:.4f} µg/kg\n")


    print("--- Step 2: Calculate PFHxS Concentration in Soil Solution (Cw) ---")
    kd_L_kg = koc * f_oc
    print(f"Soil-water partition coefficient (Kd) = Koc * f_oc = {koc:.2f} L/kg * {f_oc} = {kd_L_kg:.4f} L/kg")

    # Partitioning equation: Cw = (Cs * rho_b) / (theta_w + Kd * rho_b)
    cw_numerator = cs_ug_kg * rho_b_kg_L
    cw_denominator = theta_w + (kd_L_kg * rho_b_kg_L)
    cw_ug_L = cw_numerator / cw_denominator
    print(f"Concentration in soil solution (Cw) = ({cs_ug_kg:.4f} * {rho_b_kg_L}) / ({theta_w} + {kd_L_kg:.4f} * {rho_b_kg_L}) = {cw_ug_L:.6f} µg/L\n")


    print("--- Step 3: Calculate Total Daily Intake (DI) ---")
    # Fruit Intake
    c_fruit_ug_kg = cw_ug_L * tscf * puf_fruit
    di_fruit_ug_day = c_fruit_ug_kg * ir_fruit_kg_day * bf_fruit
    print("Daily Intake from Fruits (DI_fruit) = (Cw * TSCF * PUF_fruit) * IR_fruit * BF_fruit")
    print(f"DI_fruit = ({cw_ug_L:.6f} µg/L * {tscf} * {puf_fruit}) * {ir_fruit_kg_day} kg/day * {bf_fruit}")
    print(f"DI_fruit = {c_fruit_ug_kg:.6f} µg/kg * {ir_fruit_kg_day} kg/day * {bf_fruit} = {di_fruit_ug_day:.6f} µg/day\n")

    # Legume Intake
    c_legume_ug_kg = cw_ug_L * tscf * puf_legume
    di_legume_ug_day = c_legume_ug_kg * ir_legume_kg_day * bf_legume
    print("Daily Intake from Legumes (DI_legume) = (Cw * TSCF * PUF_legume) * IR_legume * BF_legume")
    print(f"DI_legume = ({cw_ug_L:.6f} µg/L * {tscf} * {puf_legume}) * {ir_legume_kg_day} kg/day * {bf_legume}")
    print(f"DI_legume = {c_legume_ug_kg:.6f} µg/kg * {ir_legume_kg_day} kg/day * {bf_legume} = {di_legume_ug_day:.6f} µg/day\n")
    
    # Total Intake
    di_total_ug_day = di_fruit_ug_day + di_legume_ug_day
    print(f"Total Daily Intake (DI_total) = {di_fruit_ug_day:.6f} µg/day + {di_legume_ug_day:.6f} µg/day = {di_total_ug_day:.6f} µg/day\n")


    print("--- Step 4: Calculate Hazard Quotient (HQ) ---")
    # Calculate Dose
    dose_ug_kg_day = di_total_ug_day / body_weight_kg
    print("Dose = Total Daily Intake / Body Weight")
    print(f"Dose = {di_total_ug_day:.6f} µg/day / {body_weight_kg} kg = {dose_ug_kg_day:.8f} µg/kg/day\n")

    # Calculate HQ
    hq = dose_ug_kg_day / rfd_ug_kg_day
    print("Hazard Quotient (HQ) = Dose / Reference Dose (RfD)")
    print(f"HQ = {dose_ug_kg_day:.8f} µg/kg/day / {rfd_ug_kg_day} µg/kg/day")
    print(f"Final Hazard Quotient = {hq:.6f}")
    
    return hq

if __name__ == '__main__':
    final_hq = calculate_hazard_quotient()
