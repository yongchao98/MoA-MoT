import sys

# Step 1: Define all the given variables from the problem description.
# --- Soil and Contaminant Properties ---
area_m2 = 250000
depth_m = 0.6
bulk_density_kg_m3 = 1500
foam_volume_L = 1000
pfhxs_concentration_in_foam_ug_L = 1000000

# --- Exposure and Consumption Details ---
body_weight_kg = 80
reference_dose_ug_kg_day = 0.02

# Fruit parameters
fruit_intake_g_day = 300
fruit_bioavailability = 0.5
fruit_plant_uptake_factor = 0.1
fruit_tscf = 5

# Legume parameters
legume_intake_g_day = 50
legume_bioavailability = 0.3
legume_plant_uptake_factor = 0.2
legume_tscf = 5

# Step 2: Calculate the total mass of PFHxS applied and the total mass of the soil.
total_pfhxs_mass_ug = foam_volume_L * pfhxs_concentration_in_foam_ug_L
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * bulk_density_kg_m3

# Step 3: Calculate the concentration of PFHxS in the soil (C_soil).
c_soil_ug_kg = total_pfhxs_mass_ug / soil_mass_kg

# Step 4: Calculate the PFHxS concentration in fruits and legumes.
# Equation: C_produce = C_soil * Plant_Uptake_Factor * Transpiration_Stream_Concentration_Factor
c_fruit_ug_kg = c_soil_ug_kg * fruit_plant_uptake_factor * fruit_tscf
c_legume_ug_kg = c_soil_ug_kg * legume_plant_uptake_factor * legume_tscf

# Step 5: Calculate the daily intake of PFHxS from each food source.
# Convert intake rate from g/day to kg/day.
fruit_intake_kg_day = fruit_intake_g_day / 1000
legume_intake_kg_day = legume_intake_g_day / 1000

# Equation: Daily_Intake = C_produce * Intake_Rate * Bioavailability
daily_intake_fruit_ug_day = c_fruit_ug_kg * fruit_intake_kg_day * fruit_bioavailability
daily_intake_legume_ug_day = c_legume_ug_kg * legume_intake_kg_day * legume_bioavailability

# Step 6: Calculate the total chronic daily intake (CDI) normalized by body weight.
total_daily_intake_ug_day = daily_intake_fruit_ug_day + daily_intake_legume_ug_day
cdi_ug_kg_day = total_daily_intake_ug_day / body_weight_kg

# Step 7: Calculate the Hazard Quotient (HQ).
# Equation: HQ = CDI / Reference_Dose
hazard_quotient = cdi_ug_kg_day / reference_dose_ug_kg_day

# Step 8: Print the final output showing the numbers in the final equation.
print(f"The Chronic Daily Intake (CDI) is {cdi_ug_kg_day:.6f} µg/kg/day.")
print(f"The Reference Dose (RfD) is {reference_dose_ug_kg_day} µg/kg/day.")
print("\nHazard Quotient Equation:")
print(f"Hazard Quotient = CDI / Reference Dose")
print(f"Hazard Quotient = {cdi_ug_kg_day} / {reference_dose_ug_kg_day} = {hazard_quotient}")

# Use this to output the final answer for the platform.
# print(f'<<<{hazard_quotient}>>>', file=sys.stderr)