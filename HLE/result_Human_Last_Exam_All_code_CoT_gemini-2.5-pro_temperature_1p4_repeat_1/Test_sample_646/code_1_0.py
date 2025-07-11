import math

# --- 1. Define Given Variables ---
# Contaminant and Soil Properties
foam_volume_L = 1000  # L
pfhxs_conc_in_foam_ug_per_L = 1000000  # μg/L
area_m2 = 250000  # m²
depth_m = 0.6  # m
soil_foc = 0.03  # 3% organic carbon content
soil_theta_w = 0.35  # volumetric water content (L water/L soil)
soil_rho_b_kg_m3 = 1500  # kg/m³

# Human and Food Exposure
body_weight_kg = 80  # kg

# Fruits
fruit_intake_g_day = 300  # g/day
fruit_bf = 0.5  # bioavailability factor
fruit_puf = 0.1  # plant uptake factor
fruit_tscf = 5  # transpiration stream concentration factor

# Legumes
legume_intake_g_day = 50  # g/day
legume_bf = 0.3  # bioavailability factor
legume_puf = 0.2  # plant uptake factor
legume_tscf = 5  # transpiration stream concentration factor

# Toxicological Data
rfd_ug_kg_day = 0.02  # Reference Dose in μg/kg/day

# --- 2. Calculations ---

# Step 1: Calculate total mass of PFHxS applied
total_pfhxs_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L

# Step 2: Calculate total mass of affected soil
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * soil_rho_b_kg_m3

# Step 3: Calculate PFHxS concentration in soil (C_soil)
c_soil_ug_kg = total_pfhxs_ug / soil_mass_kg

# Step 4: Calculate PFHxS concentration in soil solution (C_solution)
# The Koc for PFHxS is not given, a literature value of log(Koc)=2.8 is assumed.
koc_L_kg = 10**2.8
kd_L_kg = koc_L_kg * soil_foc

# Convert soil bulk density from kg/m³ to kg/L for the equation
soil_rho_b_kg_L = soil_rho_b_kg_m3 / 1000

# Calculate C_solution using the partitioning equation
c_solution_ug_L = c_soil_ug_kg / (kd_L_kg + (soil_theta_w / soil_rho_b_kg_L))

# Step 5: Calculate PFHxS concentration in foods (C_food)
# We assume C_food = C_solution * TSCF * PUF, and that 1L of produce weighs 1kg.
c_fruit_ug_kg = c_solution_ug_L * fruit_tscf * fruit_puf
c_legume_ug_kg = c_solution_ug_L * legume_tscf * legume_puf

# Step 6: Calculate Chronic Daily Intake (CDI)
# Convert food intake from g/day to kg/day
fruit_intake_kg_day = fruit_intake_g_day / 1000
legume_intake_kg_day = legume_intake_g_day / 1000

# Calculate absorbed dose from each food source
absorbed_dose_fruit_ug_day = c_fruit_ug_kg * fruit_intake_kg_day * fruit_bf
absorbed_dose_legume_ug_day = c_legume_ug_kg * legume_intake_kg_day * legume_bf

# Calculate total absorbed dose and CDI
total_absorbed_dose_ug_day = absorbed_dose_fruit_ug_day + absorbed_dose_legume_ug_day
cdi_ug_kg_day = total_absorbed_dose_ug_day / body_weight_kg

# Step 7: Calculate the Hazard Quotient (HQ)
hazard_quotient = cdi_ug_kg_day / rfd_ug_kg_day

# --- 3. Output the Final Calculation ---
print("This calculation determines the Hazard Quotient (HQ) based on the provided data.")
print("The final calculation is HQ = Chronic Daily Intake (CDI) / Reference Dose (RfD).\n")

print("Chronic Daily Intake (CDI) Calculation:")
print(f"  Intake from Fruits = C_fruit ({c_fruit_ug_kg:.4f} µg/kg) * IR_fruit ({fruit_intake_kg_day} kg/day) * BF_fruit ({fruit_bf}) = {absorbed_dose_fruit_ug_day:.4f} µg/day")
print(f"  Intake from Legumes = C_legume ({c_legume_ug_kg:.4f} µg/kg) * IR_legume ({legume_intake_kg_day} kg/day) * BF_legume ({legume_bf}) = {absorbed_dose_legume_ug_day:.4f} µg/day")
print(f"  Total Daily Intake = {total_absorbed_dose_ug_day:.4f} µg/day")
print(f"  CDI = Total Intake ({total_absorbed_dose_ug_day:.4f} µg/day) / Body Weight ({body_weight_kg} kg) = {cdi_ug_kg_day:.6f} µg/kg/day\n")

print("Hazard Quotient (HQ) Calculation:")
print(f"Hazard Quotient = CDI ({cdi_ug_kg_day:.6f} µg/kg/day) / RfD ({rfd_ug_kg_day} µg/kg/day)")
print(f"Hazard Quotient = {hazard_quotient:.4f}")
