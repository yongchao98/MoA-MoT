import sys

# Step 1: Define initial parameters from the problem description
# --- Soil and Contamination Data ---
soil_area_m2 = 250000
soil_depth_m = 0.6
soil_bulk_density_kg_m3 = 1500
foam_volume_L = 1000
pfhxs_concentration_in_foam_ug_L = 1000000

# --- Exposure Data ---
body_weight_kg = 80
reference_dose_ug_kg_day = 0.02

# --- Fruit Consumption Data ---
fruit_intake_g_day = 300
fruit_puf = 0.1 # Plant Uptake Factor for fruit

# --- Legume Consumption Data ---
legume_intake_g_day = 50
legume_puf = 0.2 # Plant Uptake Factor for legumes

# Note: Other parameters like organic carbon, water content, TSCF, and bioavailability
# are considered distractors as the most direct calculation path uses the Plant Uptake Factor (PUF).
# Bioavailability is typically not included when calculating a Hazard Quotient, as the
# Reference Dose (RfD) is based on external exposure (intake).

print("Step 1: Calculate the concentration of PFHxS in the soil (C_soil).\n")

# Calculate total mass of PFHxS applied
total_pfhxs_mass_ug = foam_volume_L * pfhxs_concentration_in_foam_ug_L
print(f"Total mass of PFHxS applied = {foam_volume_L} L * {pfhxs_concentration_in_foam_ug_L} ug/L = {total_pfhxs_mass_ug:,.0f} ug")

# Calculate total mass of contaminated soil
soil_volume_m3 = soil_area_m2 * soil_depth_m
soil_mass_kg = soil_volume_m3 * soil_bulk_density_kg_m3
print(f"Total mass of soil = ({soil_area_m2} m^2 * {soil_depth_m} m) * {soil_bulk_density_kg_m3} kg/m^3 = {soil_mass_kg:,.0f} kg")

# Calculate concentration in soil
c_soil_ug_kg = total_pfhxs_mass_ug / soil_mass_kg
print(f"Soil Concentration (C_soil) = {total_pfhxs_mass_ug:,.0f} ug / {soil_mass_kg:,.0f} kg = {c_soil_ug_kg:.4f} ug/kg")

print("\n" + "="*50 + "\n")

print("Step 2: Calculate the concentration of PFHxS in the consumed produce.\n")

# Calculate concentration in fruit
c_fruit_ug_kg = c_soil_ug_kg * fruit_puf
print(f"Fruit Concentration = C_soil * PUF_fruit = {c_soil_ug_kg:.4f} ug/kg * {fruit_puf} = {c_fruit_ug_kg:.4f} ug/kg")

# Calculate concentration in legumes
c_legume_ug_kg = c_soil_ug_kg * legume_puf
print(f"Legume Concentration = C_soil * PUF_legume = {c_soil_ug_kg:.4f} ug/kg * {legume_puf} = {c_legume_ug_kg:.4f} ug/kg")

print("\n" + "="*50 + "\n")

print("Step 3: Calculate the Total Daily Intake (TDI) of PFHxS.\n")

# Convert intake rates from g/day to kg/day
fruit_intake_kg_day = fruit_intake_g_day / 1000
legume_intake_kg_day = legume_intake_g_day / 1000

# Calculate Average Daily Dose (ADD) from fruit
add_fruit = (c_fruit_ug_kg * fruit_intake_kg_day) / body_weight_kg
print(f"Daily Dose from Fruit = ({c_fruit_ug_kg:.4f} ug/kg * {fruit_intake_kg_day} kg/day) / {body_weight_kg} kg = {add_fruit:.6f} ug/kg/day")

# Calculate Average Daily Dose (ADD) from legumes
add_legume = (c_legume_ug_kg * legume_intake_kg_day) / body_weight_kg
print(f"Daily Dose from Legumes = ({c_legume_ug_kg:.4f} ug/kg * {legume_intake_kg_day} kg/day) / {body_weight_kg} kg = {add_legume:.6f} ug/kg/day")

# Calculate Total Daily Intake (TDI)
total_daily_intake = add_fruit + add_legume
print(f"\nTotal Daily Intake (TDI) = {add_fruit:.6f} + {add_legume:.6f} = {total_daily_intake:.6f} ug/kg/day")

print("\n" + "="*50 + "\n")

print("Step 4: Calculate the Hazard Quotient (HQ).\n")

# Calculate Hazard Quotient
hazard_quotient = total_daily_intake / reference_dose_ug_kg_day
print("The final equation is: HQ = TDI / RfD")
print(f"Hazard Quotient = {total_daily_intake:.6f} ug/kg/day / {reference_dose_ug_kg_day} ug/kg/day")

# Print the final result
print(f"\nFinal Hazard Quotient (HQ) = {hazard_quotient:.4f}")

# Final Answer Block
# Redirecting print to stderr to not interfere with final answer format
original_stdout = sys.stdout
sys.stdout = sys.stderr
print(f"\n<<<The final answer is {hazard_quotient:.4f}>>>")
sys.stdout = original_stdout
<<<0.1111>>>