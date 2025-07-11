import math

# Step 1: Define initial parameters from the problem description
# Contamination Data
foam_volume_L = 1000  # L
pfhxs_conc_in_foam_ug_per_L = 1000000  # μg/L

# Site Data
area_m2 = 250000  # m²
depth_m = 0.6  # m
bulk_density_kg_per_m3 = 1500  # kg/m³
organic_carbon_fraction_f_oc = 0.03
volumetric_water_content_theta_w = 0.35

# Exposure Data
body_weight_kg = 80  # kg

# Food Consumption Data - Fruits
fruit_intake_g_per_day = 300  # g/day
fruit_baf = 0.5  # Bioavailability factor
fruit_tscf = 5   # Transpiration stream concentration factor

# Food Consumption Data - Legumes
legume_intake_g_per_day = 50  # g/day
legume_baf = 0.3 # Bioavailability factor
legume_tscf = 5  # Transpiration stream concentration factor

# Toxicological Data
reference_dose_rfd_ug_per_kg_day = 0.02 # μg/kg/day

# Assumption: The organic carbon-water partition coefficient (Koc) for PFHxS is not provided.
# A typical literature value for log Koc is 2.8.
# Koc = 10^2.8 L/kg
koc_L_per_kg = 10**2.8

print("--- Hazard Quotient Calculation ---")

# Step 2: Calculate total mass of PFHxS and soil, and concentration in soil (Cs)
print("\n[Step 1: Calculating Soil Concentration (Cs)]")
total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L
print(f"Total mass of PFHxS applied: {total_pfhxs_mass_ug:,.0f} μg")

soil_volume_m3 = area_m2 * depth_m
print(f"Volume of contaminated soil: {soil_volume_m3:,.0f} m³")

soil_mass_kg = soil_volume_m3 * bulk_density_kg_per_m3
print(f"Mass of contaminated soil: {soil_mass_kg:,.0f} kg")

cs_ug_per_kg = total_pfhxs_mass_ug / soil_mass_kg
print(f"Concentration in soil (Cs): {cs_ug_per_kg:.4f} μg/kg")

# Step 3: Calculate concentration in soil solution (Cw)
print("\n[Step 2: Calculating Soil Solution Concentration (Cw)]")
# Convert bulk density to kg/L for use in the partition equation (1 m³ = 1000 L)
bulk_density_kg_per_L = bulk_density_kg_per_m3 / 1000

# Calculate soil-water partition coefficient (Kd)
kd_L_per_kg = koc_L_per_kg * organic_carbon_fraction_f_oc
print(f"Assumed Koc value: {koc_L_per_kg:.2f} L/kg")
print(f"Calculated soil-water partition coefficient (Kd): {kd_L_per_kg:.2f} L/kg")

# Using the formula: Cw = (Cs * ρb) / (θw + Kd * ρb)
numerator_cw = cs_ug_per_kg * bulk_density_kg_per_L
denominator_cw = volumetric_water_content_theta_w + (kd_L_per_kg * bulk_density_kg_per_L)
cw_ug_per_L = numerator_cw / denominator_cw
print(f"Concentration in soil solution (Cw): {cw_ug_per_L:.4f} μg/L")

# Step 4: Calculate concentration in plants (C_plant)
print("\n[Step 3: Calculating Plant Concentration (C_plant)]")
# C_plant = Cw * TSCF. Assuming 1 L of plant tissue is ~1 kg.
c_plant_fruit_ug_per_kg = cw_ug_per_L * fruit_tscf
c_plant_legume_ug_per_kg = cw_ug_per_L * legume_tscf
print(f"Concentration in fruits (C_plant_fruit): {c_plant_fruit_ug_per_kg:.4f} μg/kg")
print(f"Concentration in legumes (C_plant_legume): {c_plant_legume_ug_per_kg:.4f} μg/kg")

# Step 5: Calculate Chronic Daily Intake (CDI)
print("\n[Step 4: Calculating Chronic Daily Intake (CDI)]")
# Convert intake rates from g/day to kg/day
fruit_intake_kg_per_day = fruit_intake_g_per_day / 1000
legume_intake_kg_per_day = legume_intake_g_per_day / 1000

# CDI = (C_plant * IR * BAF) / BW
cdi_fruit = (c_plant_fruit_ug_per_kg * fruit_intake_kg_per_day * fruit_baf) / body_weight_kg
cdi_legume = (c_plant_legume_ug_per_kg * legume_intake_kg_per_day * legume_baf) / body_weight_kg
print(f"Chronic Daily Intake from fruits (CDI_fruit): {cdi_fruit:.6f} μg/kg/day")
print(f"Chronic Daily Intake from legumes (CDI_legume): {cdi_legume:.6f} μg/kg/day")

total_cdi = cdi_fruit + cdi_legume
print(f"Total Chronic Daily Intake (Total CDI): {total_cdi:.6f} μg/kg/day")

# Step 6: Calculate Hazard Quotient (HQ)
print("\n[Step 5: Calculating the Final Hazard Quotient (HQ)]")
# HQ = Total CDI / RfD
hq = total_cdi / reference_dose_rfd_ug_per_kg_day

print(f"Hazard Quotient (HQ) = Total CDI / Reference Dose (RfD)")
print(f"HQ = {total_cdi:.6f} μg/kg/day / {reference_dose_rfd_ug_per_kg_day} μg/kg/day")
print(f"Final Hazard Quotient (HQ): {hq:.4f}")
print("<<<" + str(round(hq, 4)) + ">>>")