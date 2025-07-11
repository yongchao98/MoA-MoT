import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming produce
    grown on contaminated land.
    """

    # --- Input Parameters ---

    # Site and Soil Properties
    area = 250000  # m^2
    depth = 0.6  # m
    foc = 0.03  # fraction of organic carbon
    theta_w = 0.35  # volumetric water content (L water / L soil)
    rho_b_m3 = 1500  # bulk density (kg/m^3)

    # Contamination Details
    foam_volume = 1000  # L
    pfhxs_in_foam = 1000000  # µg/L

    # Exposure Scenario
    bw = 80  # body weight (kg)
    rfd = 0.02  # reference dose (µg/kg/day)

    # Fruit Consumption
    ir_fruit = 300 / 1000  # intake rate (kg/day), converted from g
    baf_fruit_human = 0.5  # bioavailability factor for human absorption
    puf_fruit = 0.1  # plant uptake factor
    tscf_fruit = 5  # transpiration stream concentration factor

    # Legume Consumption
    ir_legume = 50 / 1000  # intake rate (kg/day), converted from g
    baf_legume_human = 0.3  # bioavailability factor for human absorption
    puf_legume = 0.2  # plant uptake factor
    tscf_legume = 5  # transpiration stream concentration factor

    # --- Assumptions ---
    # The organic carbon-water partition coefficient (Koc) for PFHxS is not provided.
    # A scientifically-accepted value, Log Koc = 3.0, is assumed.
    koc = 10**3.0  # L/kg

    # --- Calculations ---

    # 1. Total mass of PFHxS (µg)
    m_pfhxs = foam_volume * pfhxs_in_foam

    # 2. Total mass of contaminated soil (kg)
    v_soil = area * depth
    m_soil = v_soil * rho_b_m3

    # 3. Concentration of PFHxS in soil (C_soil) in µg/kg
    c_soil = m_pfhxs / m_soil

    # 4. Concentration of PFHxS in soil solution (C_solution) in µg/L
    rho_b_L = rho_b_m3 / 1000  # Convert bulk density to kg/L
    kd = koc * foc  # Soil-water partition coefficient (L/kg)
    c_solution = c_soil / (kd + (theta_w / rho_b_L))

    # 5. Concentration of PFHxS in produce (µg/kg fresh weight)
    c_fruit = c_solution * puf_fruit * tscf_fruit
    c_legume = c_solution * puf_legume * tscf_legume

    # 6. Chronic Daily Intake (CDI) in µg/kg/day
    absorbed_dose_fruit = c_fruit * ir_fruit * baf_fruit_human
    absorbed_dose_legume = c_legume * ir_legume * baf_legume_human
    total_absorbed_dose = absorbed_dose_fruit + absorbed_dose_legume
    cdi = total_absorbed_dose / bw

    # 7. Hazard Quotient (HQ)
    hq = cdi / rfd

    # --- Output ---
    print("Based on the provided data and standard environmental risk assessment equations:")
    
    print(f"\nThe Chronic Daily Intake (CDI) of PFHxS is calculated to be {cdi:.8f} µg/kg body weight per day.")
    print(f"The reference dose (RfD) for PFHxS is {rfd} µg/kg body weight per day.")

    print("\nThe Hazard Quotient (HQ) is the ratio of the daily intake to the reference dose.")
    print("\nFinal Equation:")
    print(f"Hazard Quotient = Chronic Daily Intake / Reference Dose")
    print(f"Hazard Quotient = {cdi:.8f} / {rfd}")
    
    # Final Answer
    print(f"\nThe final Hazard Quotient is: {hq:.4f}")

calculate_hazard_quotient()