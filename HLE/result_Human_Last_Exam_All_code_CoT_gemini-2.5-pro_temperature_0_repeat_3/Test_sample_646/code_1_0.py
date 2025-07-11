import math

# --- Given Data ---

# PFHxS Contamination
volume_foam_L = 1000  # L
conc_in_foam_ug_per_L = 1000000  # μg/L

# Soil Properties
area_soil_m2 = 250000  # m²
depth_soil_m = 0.6  # m
organic_carbon_fraction = 0.03  # 3%
volumetric_water_content = 0.35  # L water/L soil
bulk_density_kg_per_m3 = 1500  # kg/m³

# Exposure Scenario (Man)
body_weight_kg = 80  # kg

# Fruit Consumption
intake_rate_fruit_g_day = 300  # g/day
baf_fruit = 0.5
puf_fruit = 0.1
tscf_fruit = 5

# Legume Consumption
intake_rate_legume_g_day = 50  # g/day
baf_legume = 0.3
puf_legume = 0.2
tscf_legume = 5

# Toxicological Data
reference_dose_ug_kg_day = 0.02  # μg/kg body weight per day

# --- Calculations ---

# Step 1: Calculate total mass of PFHxS applied
total_pfhxs_mass_ug = volume_foam_L * conc_in_foam_ug_per_L

# Step 2: Calculate total mass of contaminated soil
volume_soil_m3 = area_soil_m2 * depth_soil_m
total_soil_mass_kg = volume_soil_m3 * bulk_density_kg_per_m3

# Step 3: Calculate concentration of PFHxS in soil (C_soil)
conc_soil_ug_per_kg = total_pfhxs_mass_ug / total_soil_mass_kg

# Step 4: Calculate concentration of PFHxS in soil solution (C_solution)
# Assumption: Using a literature value for the organic carbon partition coefficient (K_oc) for PFHxS.
# log K_oc = 2.8, so K_oc = 10^2.8
k_oc_L_per_kg = 10**2.8
# Soil-water partition coefficient (K_d)
k_d_L_per_kg = k_oc_L_per_kg * organic_carbon_fraction
# Convert bulk density from kg/m³ to kg/L for the formula (1 m³ = 1000 L)
bulk_density_kg_per_L = bulk_density_kg_per_m3 / 1000
# Calculate C_solution
conc_solution_ug_per_L = conc_soil_ug_per_kg / (k_d_L_per_kg + (volumetric_water_content / bulk_density_kg_per_L))

# Step 5: Calculate concentration of PFHxS in plants (C_plant)
# The product of BAF, PUF, and TSCF is assumed to have units of L/kg to yield C_plant in ug/kg
conc_fruit_ug_per_kg = conc_solution_ug_per_L * baf_fruit * puf_fruit * tscf_fruit
conc_legume_ug_per_kg = conc_solution_ug_per_L * baf_legume * puf_legume * tscf_legume

# Step 6: Calculate daily intake (DI)
# Convert intake rates from g/day to kg/day
intake_rate_fruit_kg_day = intake_rate_fruit_g_day / 1000
intake_rate_legume_kg_day = intake_rate_legume_g_day / 1000

# Daily intake from fruits and legumes
daily_intake_fruit_ug_kg_day = (conc_fruit_ug_per_kg * intake_rate_fruit_kg_day) / body_weight_kg
daily_intake_legume_ug_kg_day = (conc_legume_ug_per_kg * intake_rate_legume_kg_day) / body_weight_kg

# Total Daily Intake (TDI)
total_daily_intake_ug_kg_day = daily_intake_fruit_ug_kg_day + daily_intake_legume_ug_kg_day

# Step 7: Calculate Hazard Quotient (HQ)
hazard_quotient = total_daily_intake_ug_kg_day / reference_dose_ug_kg_day

# --- Final Output ---
print("This calculation assumes a log K_oc of 2.8 for PFHxS, a common literature value.\n")
print("Calculation of the Hazard Quotient (HQ):")
print("HQ = (Total Daily Intake) / (Reference Dose)")
print("HQ = (Daily Intake from Fruits + Daily Intake from Legumes) / (Reference Dose)")

# Constructing the final equation string with all the numbers
# Daily Intake = (C_plant * IR) / BW = (C_solution * BAF * PUF * TSCF * IR) / BW
fruit_intake_str = f"(({conc_solution_ug_per_L:.4f} * {baf_fruit} * {puf_fruit} * {tscf_fruit} * {intake_rate_fruit_kg_day}) / {body_weight_kg})"
legume_intake_str = f"(({conc_solution_ug_per_L:.4f} * {baf_legume} * {puf_legume} * {tscf_legume} * {intake_rate_legume_kg_day}) / {body_weight_kg})"

print(f"\nHazard Quotient = [ {fruit_intake_str} + {legume_intake_str} ] / {reference_dose_ug_kg_day}")
print(f"Hazard Quotient = [ {daily_intake_fruit_ug_kg_day:.8f} + {daily_intake_legume_ug_kg_day:.8f} ] / {reference_dose_ug_kg_day}")
print(f"Hazard Quotient = {total_daily_intake_ug_kg_day:.8f} / {reference_dose_ug_kg_day}")
print(f"Hazard Quotient = {hazard_quotient:.4f}")
