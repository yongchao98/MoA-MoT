import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming contaminated produce.
    """
    # --- Given Data ---
    # Contamination Event
    area = 250000  # m²
    soil_depth = 0.6  # m
    foam_volume = 1000  # L
    foam_concentration = 1000000  # µg/L

    # Soil Properties
    foc_percent = 3.0  # %
    theta_w = 0.35  # L water / L soil
    rho_b = 1500  # kg/m³

    # Exposure Scenario
    body_weight = 80  # kg
    ref_dose = 0.02  # µg/kg body weight per day

    # Fruit Consumption
    intake_rate_fruit = 300 / 1000  # kg/day (converted from g)
    bioavailability_fruit = 0.5
    tscf_fruit = 5

    # Legume Consumption
    intake_rate_legume = 50 / 1000  # kg/day (converted from g)
    bioavailability_legume = 0.3
    tscf_legume = 5 # Assuming the same TSCF as it's the same chemical and conditions

    # --- Step 1: Calculate Total Mass of PFHxS ---
    total_pfhxs_mass = foam_volume * foam_concentration
    print("Step 1: Calculate Total Mass of PFHxS")
    print(f"Total PFHxS Mass = {foam_volume} L * {foam_concentration} µg/L = {total_pfhxs_mass:,.0f} µg\n")

    # --- Step 2: Calculate Total Mass of Soil ---
    soil_volume_m3 = area * soil_depth
    total_soil_mass_kg = soil_volume_m3 * rho_b
    print("Step 2: Calculate Total Mass of Soil")
    print(f"Total Soil Mass = ({area} m² * {soil_depth} m) * {rho_b} kg/m³ = {total_soil_mass_kg:,.0f} kg\n")

    # --- Step 3: Calculate PFHxS Concentration in Soil (C_soil) ---
    c_soil = total_pfhxs_mass / total_soil_mass_kg
    print("Step 3: Calculate PFHxS Concentration in Soil (C_soil)")
    print(f"C_soil = {total_pfhxs_mass:,.0f} µg / {total_soil_mass_kg:,.0f} kg = {c_soil:.4f} µg/kg\n")

    # --- Step 4: Calculate PFHxS Concentration in Soil Pore Water (C_pw) ---
    # Koc for PFHxS is not given, a literature value of 10^2.8 L/kg is used.
    log_koc = 2.8
    koc = math.pow(10, log_koc)
    foc = foc_percent / 100
    kd = foc * koc
    
    # Convert bulk density to kg/L (1 m³ = 1000 L)
    rho_b_kg_per_L = rho_b / 1000

    # C_pw = C_soil / (Kd + θw / ρb)
    c_pw = c_soil / (kd + (theta_w / rho_b_kg_per_L))
    print("Step 4: Calculate PFHxS Concentration in Soil Pore Water (C_pw)")
    print(f"Kd (Soil-Water Partition Coefficient) = {foc} * {koc:.2f} L/kg = {kd:.2f} L/kg")
    print(f"C_pw = {c_soil:.4f} µg/kg / ({kd:.2f} L/kg + {theta_w} / {rho_b_kg_per_L} kg/L) = {c_pw:.4f} µg/L\n")

    # --- Step 5: Calculate PFHxS Concentration in Produce (C_plant) ---
    # Assuming density of produce is 1 kg/L, so C_plant in µg/L is equivalent to µg/kg.
    # The 'plant uptake factor' is considered redundant as TSCF is a more specific parameter for this calculation.
    c_plant = c_pw * tscf_fruit # TSCF is the same for both
    print("Step 5: Calculate PFHxS Concentration in Produce (C_plant)")
    print(f"C_plant = C_pw * TSCF = {c_pw:.4f} µg/L * {tscf_fruit} = {c_plant:.4f} µg/kg\n")

    # --- Step 6: Calculate Estimated Daily Intake (EDI) ---
    absorbed_intake_fruit = c_plant * intake_rate_fruit * bioavailability_fruit
    absorbed_intake_legume = c_plant * intake_rate_legume * bioavailability_legume
    total_absorbed_intake = absorbed_intake_fruit + absorbed_intake_legume
    edi = total_absorbed_intake / body_weight
    print("Step 6: Calculate Estimated Daily Intake (EDI)")
    print(f"EDI = ( (C_plant * IR_fruit * AF_fruit) + (C_plant * IR_legume * AF_legume) ) / Body Weight")
    print(f"EDI = ( ({c_plant:.4f} * {intake_rate_fruit} * {bioavailability_fruit}) + ({c_plant:.4f} * {intake_rate_legume} * {bioavailability_legume}) ) / {body_weight}")
    print(f"EDI = {edi:.6f} µg/kg/day\n")

    # --- Step 7: Calculate Hazard Quotient (HQ) ---
    hq = edi / ref_dose
    print("Step 7: Calculate the Hazard Quotient (HQ)")
    print(f"Hazard Quotient = EDI / Reference Dose")
    print(f"Hazard Quotient = {edi:.6f} µg/kg/day / {ref_dose} µg/kg/day = {hq:.4f}")

    return hq

if __name__ == '__main__':
    hazard_quotient = calculate_hazard_quotient()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<< {hazard_quotient:.4f} >>>")