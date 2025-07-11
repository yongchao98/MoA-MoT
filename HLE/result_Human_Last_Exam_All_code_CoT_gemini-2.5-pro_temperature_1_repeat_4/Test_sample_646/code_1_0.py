import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming
    contaminated produce based on the provided scenario.
    """
    # Step 0: Define constants and variables from the problem
    # Soil properties
    area_m2 = 250000
    depth_m = 0.6
    f_oc = 0.03
    theta_w = 0.35
    rho_b_kg_m3 = 1500

    # Contaminant properties
    foam_volume_L = 1000
    foam_conc_ug_L = 1000000

    # Exposure scenario
    body_weight_kg = 80

    # Fruit consumption
    fruit_intake_g_day = 300
    fruit_bioavailability = 0.5
    fruit_puf = 0.1

    # Legume consumption
    legume_intake_g_day = 50
    legume_bioavailability = 0.3
    legume_puf = 0.2

    # Common factors
    tscf = 5

    # Toxicological data
    rfd_ug_kg_day = 0.02

    # Assumed literature value for PFHxS Koc
    log_koc = 2.7
    k_oc_L_kg = 10**log_koc

    print("--- Hazard Quotient Calculation for PFHxS Exposure ---\n")

    # Step 1: Calculate the total mass of PFHxS applied
    total_pfhxs_ug = foam_volume_L * foam_conc_ug_L
    print(f"Step 1: Calculate the total mass of PFHxS applied to the soil.")
    print(f"Total PFHxS (μg) = {foam_volume_L} L * {foam_conc_ug_L:,} μg/L = {total_pfhxs_ug:,.0f} μg")
    print("-" * 50)

    # Step 2: Calculate the total mass of the affected soil
    soil_volume_m3 = area_m2 * depth_m
    soil_mass_kg = soil_volume_m3 * rho_b_kg_m3
    print(f"Step 2: Calculate the total mass of the affected soil.")
    print(f"Soil Volume (m³) = {area_m2:,} m² * {depth_m} m = {soil_volume_m3:,.0f} m³")
    print(f"Soil Mass (kg) = {soil_volume_m3:,.0f} m³ * {rho_b_kg_m3} kg/m³ = {soil_mass_kg:,.0f} kg")
    print("-" * 50)

    # Step 3: Calculate the concentration of PFHxS in the soil (C_soil)
    c_soil_ug_kg = total_pfhxs_ug / soil_mass_kg
    print(f"Step 3: Calculate the concentration of PFHxS in the soil (C_soil).")
    print(f"C_soil (μg/kg) = {total_pfhxs_ug:,.0f} μg / {soil_mass_kg:,.0f} kg = {c_soil_ug_kg:.4f} μg/kg")
    print("-" * 50)

    # Step 4: Calculate the concentration of PFHxS in the soil solution (C_sw)
    rho_b_kg_L = rho_b_kg_m3 / 1000
    kd_L_kg = k_oc_L_kg * f_oc
    c_sw_ug_L = c_soil_ug_kg / (kd_L_kg + (theta_w / rho_b_kg_L))
    print(f"Step 4: Calculate the concentration of PFHxS in the soil solution (C_sw).")
    print(f"Using an assumed log(Koc) of {log_koc}, Koc = {k_oc_L_kg:.2f} L/kg.")
    print(f"Soil-water partition coefficient (Kd) = {k_oc_L_kg:.2f} L/kg * {f_oc} = {kd_L_kg:.4f} L/kg")
    print(f"Bulk density (ρb) = {rho_b_kg_m3} kg/m³ = {rho_b_kg_L} kg/L")
    print(f"C_sw (μg/L) = {c_soil_ug_kg:.4f} μg/kg / ({kd_L_kg:.4f} L/kg + ({theta_w} L/L / {rho_b_kg_L} kg/L)) = {c_sw_ug_L:.4f} μg/L")
    print("-" * 50)

    # Step 5: Calculate the concentration of PFHxS in the consumed produce
    # Assuming 1 L of plant tissue has a mass of 1 kg, C(μg/kg) = C(μg/L)
    # Convert to μg/g by dividing by 1000
    c_fruit_ug_g = (c_sw_ug_L * tscf * fruit_puf) / 1000
    c_legume_ug_g = (c_sw_ug_L * tscf * legume_puf) / 1000
    print(f"Step 5: Calculate the concentration of PFHxS in the consumed produce.")
    print(f"C_fruit (μg/g) = ({c_sw_ug_L:.4f} μg/L * {tscf} * {fruit_puf}) / 1000 = {c_fruit_ug_g:.6f} μg/g")
    print(f"C_legume (μg/g) = ({c_sw_ug_L:.4f} μg/L * {tscf} * {legume_puf}) / 1000 = {c_legume_ug_g:.6f} μg/g")
    print("-" * 50)

    # Step 6: Calculate the total absorbed daily intake (TADI)
    adi_fruit_ug_kg_day = (c_fruit_ug_g * fruit_intake_g_day * fruit_bioavailability) / body_weight_kg
    adi_legume_ug_kg_day = (c_legume_ug_g * legume_intake_g_day * legume_bioavailability) / body_weight_kg
    tadi_ug_kg_day = adi_fruit_ug_kg_day + adi_legume_ug_kg_day
    print(f"Step 6: Calculate the total absorbed daily intake (TADI) of PFHxS.")
    print(f"Absorbed Daily Intake from fruits (μg/kg/day) = ({c_fruit_ug_g:.6f} μg/g * {fruit_intake_g_day} g/day * {fruit_bioavailability}) / {body_weight_kg} kg = {adi_fruit_ug_kg_day:.6f} μg/kg/day")
    print(f"Absorbed Daily Intake from legumes (μg/kg/day) = ({c_legume_ug_g:.6f} μg/g * {legume_intake_g_day} g/day * {legume_bioavailability}) / {body_weight_kg} kg = {adi_legume_ug_kg_day:.6f} μg/kg/day")
    print(f"Total Absorbed Daily Intake (μg/kg/day) = {adi_fruit_ug_kg_day:.6f} μg/kg/day + {adi_legume_ug_kg_day:.6f} μg/kg/day = {tadi_ug_kg_day:.6f} μg/kg/day")
    print("-" * 50)

    # Step 7: Calculate the Hazard Quotient (HQ)
    hazard_quotient = tadi_ug_kg_day / rfd_ug_kg_day
    print(f"Step 7: Calculate the final Hazard Quotient (HQ).")
    print(f"HQ = Total Absorbed Daily Intake / Reference Dose")
    print(f"HQ = {tadi_ug_kg_day:.6f} μg/kg/day / {rfd_ug_kg_day} μg/kg/day = {hazard_quotient:.4f}")
    print("-" * 50)

    return hazard_quotient

if __name__ == '__main__':
    final_hq = calculate_hazard_quotient()
    print(f"\nThe final Hazard Quotient is {final_hq:.4f}.")
    print(f"<<<{final_hq}>>>")
