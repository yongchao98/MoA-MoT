import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for a man exposed to PFHxS
    from consuming produce grown on contaminated land.
    """

    # --- 1. Define Given Data and Assumptions ---
    # Chemical (PFHxS)
    v_foam = 1000  # L (Volume of firefighting foam)
    c_foam = 1000000  # µg/L (Concentration of PFHxS in foam)

    # Soil Properties
    area_soil = 250000  # m²
    depth_soil = 0.6  # m
    f_oc = 0.03  # Organic carbon content fraction (3%)
    theta_w = 0.35  # Volumetric water content (L water/L soil)
    rho_b = 1500  # Bulk density (kg/m³)

    # Exposure Scenario (30-year-old male)
    body_weight = 80  # kg

    # Consumption Data - Fruits
    ir_fruit = 300 / 1000  # kg/day (Daily intake of fruits)
    baf_fruit = 0.5  # Bioavailability factor for fruits
    puf_fruit = 0.1  # Plant uptake factor for fruits
    tscf_fruit = 5  # Transpiration stream concentration factor for fruits

    # Consumption Data - Legumes
    ir_legume = 50 / 1000  # kg/day (Daily intake of legumes)
    baf_legume = 0.3  # Bioavailability factor for legumes
    puf_legume = 0.2  # Plant uptake factor for legumes
    tscf_legume = 5  # Transpiration stream concentration factor for legumes

    # Toxicological Data
    ref_dose = 0.02  # µg/kg body weight per day (Reference Dose for PFHxS)

    # Assumption for a missing parameter
    # The Organic Carbon Partition Coefficient (Koc) for PFHxS is not given.
    # A literature value of log Koc = 2.8 is assumed.
    log_koc = 2.8
    koc = 10**log_koc  # L/kg

    # --- 2. Step-by-Step Calculations ---

    # Total mass of PFHxS applied to the soil (µg)
    total_pfhxs_mass = v_foam * c_foam

    # Total volume and mass of contaminated soil
    total_soil_volume = area_soil * depth_soil  # m³
    total_soil_mass = total_soil_volume * rho_b  # kg

    # Soil-water partition coefficient (Kd in L/kg)
    kd = koc * f_oc

    # Total volume of water in the soil (L)
    total_water_volume_litres = total_soil_volume * theta_w * 1000

    # Concentration of PFHxS in soil solution (µg/L) using a mass balance equation:
    # Total Mass = (Mass in Water) + (Mass sorbed to Soil)
    # Total Mass = C_solution * V_water + (C_solution * Kd) * M_soil
    c_solution = total_pfhxs_mass / (total_water_volume_litres + kd * total_soil_mass)

    # Concentration of PFHxS in fruits (µg/kg)
    c_fruit = c_solution * tscf_fruit * puf_fruit

    # Concentration of PFHxS in legumes (µg/kg)
    c_legume = c_solution * tscf_legume * puf_legume

    # Daily intake of PFHxS from fruits (µg/kg BW/day)
    di_fruit = (c_fruit * ir_fruit * baf_fruit) / body_weight

    # Daily intake of PFHxS from legumes (µg/kg BW/day)
    di_legume = (c_legume * ir_legume * baf_legume) / body_weight

    # Total daily intake (TDI) from both sources
    total_daily_intake = di_fruit + di_legume

    # Final Hazard Quotient (HQ)
    hazard_quotient = total_daily_intake / ref_dose

    # --- 3. Output the Final Equation ---
    print("Hazard Quotient (HQ) = (Total Daily Intake) / (Reference Dose)")
    print("Total Daily Intake = (Daily Intake from Fruits) + (Daily Intake from Legumes)")
    print("\nCalculating the final Hazard Quotient:")
    print(f"HQ = ({di_fruit:.8f} + {di_legume:.8f}) / {ref_dose:.3f}")
    print(f"HQ = {total_daily_intake:.8f} / {ref_dose:.3f}")
    print(f"HQ = {hazard_quotient:.5f}")


if __name__ == "__main__":
    calculate_hazard_quotient()
    # The final numerical answer is approximately 0.01304
    # print("<<<0.01304>>>")