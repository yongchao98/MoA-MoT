import math

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for PFHxS exposure from consuming produce
    grown on contaminated land.
    """

    # --- Step 0: Define Constants and Given Data ---

    # Site and Soil Characteristics
    area = 250000  # m^2
    depth = 0.6  # m
    f_oc = 0.03  # organic carbon fraction (3%)
    theta_w = 0.35  # volumetric water content (L water / L soil)
    rho_b_m3 = 1500  # bulk density (kg/m^3)

    # Contamination Details
    foam_volume = 1000  # L
    pfhxs_conc_in_foam = 1000000  # ug/L
    # Assumption for PFHxS soil-organic-carbon-water partition coefficient (Koc)
    # Based on literature (e.g., US EPA), Log Koc for PFHxS is around 3.0.
    koc = 10**3  # L/kg

    # Exposure Scenario
    body_weight = 80  # kg
    ref_dose = 0.02  # ug/kg body weight per day

    # Fruit Consumption Data
    ir_fruit = 300 / 1000  # kg/day (converted from g/day)
    bf_fruit = 0.5  # bioavailability factor
    puf_fruit = 0.1  # plant uptake factor
    tscf_fruit = 5  # transpiration stream concentration factor

    # Legume Consumption Data
    ir_legume = 50 / 1000  # kg/day (converted from g/day)
    bf_legume = 0.3  # bioavailability factor
    puf_legume = 0.2  # plant uptake factor
    tscf_legume = 5  # transpiration stream concentration factor

    print("--- Intermediate Calculations ---")

    # --- Step 1: Calculate total mass of PFHxS applied ---
    total_pfhxs_mass = foam_volume * pfhxs_conc_in_foam
    print(f"1. Total mass of PFHxS applied: {total_pfhxs_mass:,.0f} µg")

    # --- Step 2: Calculate mass of contaminated soil ---
    soil_volume = area * depth
    soil_mass = soil_volume * rho_b_m3
    print(f"2. Total mass of contaminated soil: {soil_mass:,.0f} kg")

    # --- Step 3: Calculate total concentration of PFHxS in soil (Ct) ---
    ct_soil = total_pfhxs_mass / soil_mass
    print(f"3. Total concentration of PFHxS in soil (Ct): {ct_soil:.4f} µg/kg")

    # --- Step 4: Calculate concentration of PFHxS in soil porewater (Cw) ---
    rho_b_l = rho_b_m3 / 1000  # Convert bulk density to kg/L
    kd = koc * f_oc # Calculate soil-water partition coefficient (Kd)
    cw = ct_soil / ((theta_w / rho_b_l) + kd)
    print(f"4. Concentration of PFHxS in soil porewater (Cw): {cw:.4f} µg/L")

    # --- Step 5: Calculate concentration of PFHxS in produce ---
    c_fruit = cw * puf_fruit * tscf_fruit
    c_legume = cw * puf_legume * tscf_legume
    print(f"5. Concentration in fruits: {c_fruit:.4f} µg/kg | Concentration in legumes: {c_legume:.4f} µg/kg")

    # --- Step 6: Calculate the total daily intake of PFHxS ---
    di_fruit = c_fruit * ir_fruit * bf_fruit
    di_legume = c_legume * ir_legume * bf_legume
    di_total = di_fruit + di_legume
    print(f"6. Total Daily Intake (DI): {di_total:.4f} µg/day")

    # --- Step 7: Calculate the Chronic Daily Intake (CDI) ---
    cdi = di_total / body_weight
    print(f"7. Chronic Daily Intake (CDI): {cdi:.6f} µg/kg/day")

    # --- Step 8: Calculate the Hazard Quotient (HQ) ---
    hq = cdi / ref_dose
    print("\n--- Final Calculation ---")
    print(f"The Hazard Quotient is calculated by dividing the Chronic Daily Intake (CDI) by the Reference Dose (RfD).")
    print(f"Hazard Quotient = {cdi:.6f} / {ref_dose} = {hq:.5f}")
    
    # Returning the value for the final answer block
    return hq

# Run the calculation and store the result
hazard_quotient_value = calculate_hazard_quotient()
# The final answer format is handled outside the function printout for clarity
final_answer = f"<<<{hazard_quotient_value:.3f}>>>"

# The final answer part is for the platform, but for a user running the script, this print is useful.
# print(f"\nFinal Answer: {final_answer}")
# To comply with instructions, the last thing in the script's output should be the answer.
print(f"{final_answer}")

if __name__ == '__main__':
    # The functions are called when the script is run directly.
    # The output is already generated above.
    pass