import math

# Step 1: Define all the given parameters

# Contamination Details
foam_volume_L = 1000  # L
pfhxs_conc_in_foam_ug_per_L = 1000000  # μg/L

# Soil Characteristics
area_m2 = 250000  # m²
depth_m = 0.6  # m
organic_carbon_fraction = 0.03  # 3%
volumetric_water_content = 0.35  # L water/L soil
bulk_density_kg_per_m3 = 1500  # kg/m³

# Exposure Scenario (Man)
body_weight_kg = 80  # kg

# Fruit Consumption
fruit_intake_g_per_day = 300  # g/day
fruit_bioavailability = 0.5
fruit_plant_uptake_factor = 0.1
fruit_transpiration_factor = 5

# Legume Consumption
legume_intake_g_per_day = 50  # g/day
legume_bioavailability = 0.3
legume_plant_uptake_factor = 0.2
legume_transpiration_factor = 5

# Toxicological Data
reference_dose_ug_per_kg_day = 0.02  # μg/kg/day

# Step 2: Perform Calculations

# Convert units for consistency
fruit_intake_kg_per_day = fruit_intake_g_per_day / 1000
legume_intake_kg_per_day = legume_intake_g_per_day / 1000
bulk_density_kg_per_L = bulk_density_kg_per_m3 / 1000

# Total mass of PFHxS applied
total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L

# Total volume and mass of contaminated soil
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * bulk_density_kg_per_m3

# Concentration of PFHxS in bulk soil (Cs)
Cs_ug_per_kg = total_pfhxs_mass_ug / soil_mass_kg

# Concentration in soil solution (Cw)
# Using a literature value for log Koc of PFHxS = 2.8
Koc_L_per_kg = 10**2.8
Kd_L_per_kg = Koc_L_per_kg * organic_carbon_fraction
Cw_ug_per_L = Cs_ug_per_kg / (Kd_L_per_kg + (volumetric_water_content / bulk_density_kg_per_L))

# Concentration in plants (C_plant)
# BAF = Plant Uptake Factor * Transpiration Stream Concentration Factor
baf_fruit = fruit_plant_uptake_factor * fruit_transpiration_factor
C_fruit_ug_per_kg = Cw_ug_per_L * baf_fruit

baf_legume = legume_plant_uptake_factor * legume_transpiration_factor
C_legume_ug_per_kg = Cw_ug_per_L * baf_legume

# Daily Intake (DI) from each food source
DI_fruit_ug_per_day = C_fruit_ug_per_kg * fruit_intake_kg_per_day * fruit_bioavailability
DI_legume_ug_per_day = C_legume_ug_per_kg * legume_intake_kg_per_day * legume_bioavailability

# Total Daily Intake (TDI) and Chronic Daily Intake (CDI)
TDI_ug_per_day = DI_fruit_ug_per_day + DI_legume_ug_per_day
CDI_ug_per_kg_day = TDI_ug_per_day / body_weight_kg

# Hazard Quotient (HQ)
HQ = CDI_ug_per_kg_day / reference_dose_ug_per_kg_day

# Step 3: Print the final equation with calculated values
print("Calculation of the Hazard Quotient (HQ)")
print("="*40)
print(f"Daily Intake from Fruits (ug/day) = {C_fruit_ug_per_kg:.4f} ug/kg * {fruit_intake_kg_per_day} kg/day * {fruit_bioavailability} = {DI_fruit_ug_per_day:.4f} ug/day")
print(f"Daily Intake from Legumes (ug/day) = {C_legume_ug_per_kg:.4f} ug/kg * {legume_intake_kg_per_day} kg/day * {legume_bioavailability} = {DI_legume_ug_per_day:.4f} ug/day")
print(f"Total Daily Intake (ug/day) = {DI_fruit_ug_per_day:.4f} + {DI_legume_ug_per_day:.4f} = {TDI_ug_per_day:.4f} ug/day")
print(f"Chronic Daily Intake (ug/kg/day) = {TDI_ug_per_day:.4f} ug/day / {body_weight_kg} kg = {CDI_ug_per_kg_day:.6f} ug/kg/day")
print("\n--- Final Hazard Quotient Equation ---")
print(f"HQ = Chronic Daily Intake / Reference Dose")
print(f"HQ = {CDI_ug_per_kg_day:.6f} / {reference_dose_ug_per_kg_day}")
print(f"HQ = {HQ:.4f}")

print("\n--- Final Answer ---")
print(f"The Hazard Quotient (HQ) is the total daily intake of PFHxS per unit of body weight, divided by the reference dose.")
print(f"HQ = (({DI_fruit_ug_per_day:.4f} + {DI_legume_ug_per_day:.4f}) ug/day / {body_weight_kg} kg) / {reference_dose_ug_per_kg_day} ug/kg/day")
print(f"HQ = ({TDI_ug_per_day:.4f} ug/day / {body_weight_kg} kg) / {reference_dose_ug_per_kg_day} ug/kg/day")
print(f"HQ = {CDI_ug_per_kg_day:.6f} ug/kg/day / {reference_dose_ug_per_kg_day} ug/kg/day")
print(f"HQ = {HQ:.4f}")