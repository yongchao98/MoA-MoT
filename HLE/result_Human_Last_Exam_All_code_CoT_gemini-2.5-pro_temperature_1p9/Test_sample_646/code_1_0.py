import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

def calculate_hazard_quotient():
    """
    Calculates the Hazard Quotient (HQ) for exposure to PFHxS from consuming contaminated produce.
    The steps are printed to explain the process.
    """

    # --- Given Data ---
    # Contamination event
    foam_volume_L = 1000
    pfhxs_in_foam_ug_per_L = 1000000

    # Soil properties
    area_m2 = 250000
    depth_m = 0.6
    bulk_density_kg_per_m3 = 1500

    # Human exposure
    body_weight_kg = 80
    reference_dose_ug_per_kg_day = 0.02

    # Fruit consumption
    intake_fruit_g_day = 300
    bioavailability_fruit = 0.5
    puf_fruit = 0.1 # Plant Uptake Factor for fruit

    # Legume consumption
    intake_legume_g_day = 50
    bioavailability_legume = 0.3
    puf_legume = 0.2 # Plant Uptake Factor for legume

    # --- Step-by-step Calculation ---

    print("Step 1: Calculate the total mass of PFHxS applied to the soil.")
    total_pfhxs_mass_ug = foam_volume_L * pfhxs_in_foam_ug_per_L
    print(f"Total PFHxS Mass = {foam_volume_L} L * {pfhxs_in_foam_ug_per_L} μg/L = {total_pfhxs_mass_ug:,.0f} μg\n")

    print("Step 2: Calculate the total mass of the contaminated soil.")
    soil_volume_m3 = area_m2 * depth_m
    soil_mass_kg = soil_volume_m3 * bulk_density_kg_per_m3
    print(f"Soil Volume = {area_m2:,.0f} m² * {depth_m} m = {soil_volume_m3:,.0f} m³")
    print(f"Soil Mass = {soil_volume_m3:,.0f} m³ * {bulk_density_kg_per_m3} kg/m³ = {soil_mass_kg:,.0f} kg\n")

    print("Step 3: Calculate the concentration of PFHxS in the bulk soil (Cs).")
    cs_ug_per_kg = total_pfhxs_mass_ug / soil_mass_kg
    print(f"Cs = Total PFHxS Mass / Soil Mass = {total_pfhxs_mass_ug:,.0f} μg / {soil_mass_kg:,.0f} kg = {cs_ug_per_kg:.4f} μg/kg\n")

    print("Step 4: Calculate the concentration of PFHxS in the consumed produce using the Plant Uptake Factor (PUF).")
    # Convert intake rates from g/day to kg/day
    intake_fruit_kg_day = intake_fruit_g_day / 1000
    intake_legume_kg_day = intake_legume_g_day / 1000

    # Concentration in fruits
    c_fruit_ug_per_kg = cs_ug_per_kg * puf_fruit
    print(f"Concentration in Fruits = Cs * PUF_fruit = {cs_ug_per_kg:.4f} μg/kg * {puf_fruit} = {c_fruit_ug_per_kg:.4f} μg/kg")
    
    # Concentration in legumes
    c_legume_ug_per_kg = cs_ug_per_kg * puf_legume
    print(f"Concentration in Legumes = Cs * PUF_legume = {cs_ug_per_kg:.4f} μg/kg * {puf_legume} = {c_legume_ug_per_kg:.4f} μg/kg\n")

    print("Step 5: Calculate the total absorbed daily intake of PFHxS.")
    absorbed_intake_fruit_ug_day = c_fruit_ug_per_kg * intake_fruit_kg_day * bioavailability_fruit
    print(f"Absorbed from Fruits = {c_fruit_ug_per_kg:.4f} μg/kg * {intake_fruit_kg_day} kg/day * {bioavailability_fruit} = {absorbed_intake_fruit_ug_day:.4f} μg/day")
    
    absorbed_intake_legume_ug_day = c_legume_ug_per_kg * intake_legume_kg_day * bioavailability_legume
    print(f"Absorbed from Legumes = {c_legume_ug_per_kg:.4f} μg/kg * {intake_legume_kg_day} kg/day * {bioavailability_legume} = {absorbed_intake_legume_ug_day:.4f} μg/day")

    total_absorbed_intake_ug_day = absorbed_intake_fruit_ug_day + absorbed_intake_legume_ug_day
    print(f"Total Absorbed Intake = {absorbed_intake_fruit_ug_day:.4f} μg/day + {absorbed_intake_legume_ug_day:.4f} μg/day = {total_absorbed_intake_ug_day:.4f} μg/day\n")
    
    print("Step 6: Calculate the dose, normalized by body weight.")
    dose_ug_per_kg_day = total_absorbed_intake_ug_day / body_weight_kg
    print(f"Dose = Total Absorbed Intake / Body Weight = {total_absorbed_intake_ug_day:.4f} μg/day / {body_weight_kg} kg = {dose_ug_per_kg_day:.6f} μg/kg/day\n")

    print("Step 7: Calculate the final Hazard Quotient (HQ).")
    hazard_quotient = dose_ug_per_kg_day / reference_dose_ug_per_kg_day
    print(f"Hazard Quotient (HQ) = Dose / Reference Dose (RfD)")
    print(f"HQ = {dose_ug_per_kg_day:.6f} μg/kg/day / {reference_dose_ug_per_kg_day} μg/kg/day = {hazard_quotient:.4f}\n")
    
    # Final answer in required format
    print(f"The final calculated Hazard Quotient is {hazard_quotient:.4f}.")

calculate_hazard_quotient()

# Restore original stdout
sys.stdout = original_stdout
# Get the buffered output
output = buffer.getvalue()

print(output)

# Extract final number for the answer block
final_hq_value = float(output.splitlines()[-2].split()[-1].replace('.', '')) / 10000

print(f'<<<{final_hq_value}>>>')