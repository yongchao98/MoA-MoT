import math

# --- Step 1: Define Initial Variables ---
# PFHxS and Exposure
foam_volume_L = 1000
pfhxs_conc_in_foam_ug_per_L = 1000000
reference_dose_ug_per_kg_day = 0.02

# Soil Properties
area_m2 = 250000
depth_m = 0.6
f_oc = 0.03  # Organic carbon content fraction
theta_w = 0.35  # Volumetric water content (L water/L soil)
bulk_density_kg_per_m3 = 1500
k_oc_L_per_kg = 1000 # Assumed literature value for Organic Carbon-Water Partition Coefficient

# Human Receptor
body_weight_kg = 80

# Fruit Consumption
fruit_intake_g_day = 300
fruit_bioavailability = 0.5
fruit_puf = 0.1
fruit_tscf = 5

# Legume Consumption
legume_intake_g_day = 50
legume_bioavailability = 0.3
legume_puf = 0.2
legume_tscf = 5

# --- Step 2: Calculate PFHxS and Soil Mass ---
total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L
soil_volume_m3 = area_m2 * depth_m
total_soil_mass_kg = soil_volume_m3 * bulk_density_kg_per_m3

# --- Step 3: Calculate Soil and Porewater Concentrations ---
cs_ug_per_kg = total_pfhxs_mass_ug / total_soil_mass_kg
kd_L_per_kg = f_oc * k_oc_L_per_kg
cw_ug_per_L = cs_ug_per_kg / (kd_L_per_kg + theta_w)

# --- Step 4: Calculate Plant Concentrations ---
c_fruit_ug_per_kg = cw_ug_per_L * fruit_puf * fruit_tscf
c_legume_ug_per_kg = cw_ug_per_L * legume_puf * legume_tscf

# --- Step 5: Calculate Daily Intake ---
# Convert food intake from g/day to kg/day
fruit_intake_kg_day = fruit_intake_g_day / 1000
legume_intake_kg_day = legume_intake_g_day / 1000

# Daily intake formula: (C_plant * IR * BF) / BW
di_fruit_ug_per_kg_day = (c_fruit_ug_per_kg * fruit_intake_kg_day * fruit_bioavailability) / body_weight_kg
di_legume_ug_per_kg_day = (c_legume_ug_per_kg * legume_intake_kg_day * legume_bioavailability) / body_weight_kg

total_daily_intake_ug_per_kg_day = di_fruit_ug_per_kg_day + di_legume_ug_per_kg_day

# --- Step 6: Calculate Hazard Quotient ---
hazard_quotient = total_daily_intake_ug_per_kg_day / reference_dose_ug_per_kg_day

# --- Step 7: Print the Final Equation and Result ---
print("Calculation of the Hazard Quotient (HQ)")
print("-" * 40)
print(f"Total Daily Intake (DI) = DI from Fruits + DI from Legumes")
print(f"DI from Fruits = (C_fruit * IR_fruit * BF_fruit) / BW")
print(f"DI from Fruits = ({c_fruit_ug_per_kg:.6f} µg/kg * {fruit_intake_kg_day} kg/day * {fruit_bioavailability}) / {body_weight_kg} kg = {di_fruit_ug_per_kg_day:.6f} µg/kg/day")
print(f"DI from Legumes = (C_legume * IR_legume * BF_legume) / BW")
print(f"DI from Legumes = ({c_legume_ug_per_kg:.6f} µg/kg * {legume_intake_kg_day} kg/day * {legume_bioavailability}) / {body_weight_kg} kg = {di_legume_ug_per_kg_day:.6f} µg/kg/day")
print("-" * 40)
print(f"Total Daily Intake = {di_fruit_ug_per_kg_day:.6f} + {di_legume_ug_per_kg_day:.6f} = {total_daily_intake_ug_per_kg_day:.6f} µg/kg/day")
print(f"Reference Dose (RfD) = {reference_dose_ug_per_kg_day} µg/kg/day")
print("-" * 40)
print("Final Equation:")
print(f"Hazard Quotient = Total Daily Intake / Reference Dose")
print(f"Hazard Quotient = {total_daily_intake_ug_per_kg_day:.6f} / {reference_dose_ug_per_kg_day} = {hazard_quotient:.6f}")
print("-" * 40)

<<<0.008235>>>