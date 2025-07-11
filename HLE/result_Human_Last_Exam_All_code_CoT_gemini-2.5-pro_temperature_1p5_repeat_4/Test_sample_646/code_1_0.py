import math

# --- Step 0: Define Initial Parameters ---
print("--- Step 0: Initial Parameters ---")
# Contaminant info
foam_volume_L = 1000
foam_concentration_ug_per_L = 1000000

# Soil info
area_m2 = 250000
depth_m = 0.6
f_oc = 0.03  # organic carbon content fraction
theta_w = 0.35 # volumetric water content (L water / L soil)
rho_b_kg_per_m3 = 1500 # bulk density

# Human receptor info
body_weight_kg = 80

# Fruit consumption info
intake_fruit_g_day = 300
baf_fruit = 0.5
puf_fruit = 0.1 # Plant Uptake Factor (not used in TSCF model)
tscf_fruit = 5

# Legume consumption info
intake_legume_g_day = 50
baf_legume = 0.3
puf_legume = 0.2 # Plant Uptake Factor (not used in TSCF model)
tscf_legume = 5

# Toxicological info
ref_dose_ug_per_kg_day = 0.02

# --- Step 1: Calculate Total Mass of PFHxS and Soil ---
print("\n--- Step 1: Calculate Contaminant Mass and Soil Mass ---")
# Total mass of PFHxS
total_pfhxs_ug = foam_volume_L * foam_concentration_ug_per_L
print(f"Total PFHxS Mass = {foam_volume_L} L * {foam_concentration_ug_per_L} ug/L = {total_pfhxs_ug:,.0f} ug")

# Volume of soil
volume_soil_m3 = area_m2 * depth_m
print(f"Volume of Soil = {area_m2:,.0f} m^2 * {depth_m} m = {volume_soil_m3:,.0f} m^3")

# Total mass of soil
total_soil_kg = volume_soil_m3 * rho_b_kg_per_m3
print(f"Total Soil Mass = {volume_soil_m3:,.0f} m^3 * {rho_b_kg_per_m3} kg/m^3 = {total_soil_kg:,.0f} kg")

# --- Step 2: Calculate PFHxS Concentration in Soil (Bulk and Solution) ---
print("\n--- Step 2: Calculate PFHxS Concentration in Soil ---")
# Concentration in bulk soil (C_soil)
c_soil_ug_per_kg = total_pfhxs_ug / total_soil_kg
print(f"Bulk Soil Concentration (C_soil) = {total_pfhxs_ug:,.0f} ug / {total_soil_kg:,.0f} kg = {c_soil_ug_per_kg:.4f} ug/kg")

# This calculation requires a K_oc value (organic carbon-water partition coefficient) for PFHxS.
# This was not provided, so we assume a common literature value.
# log K_oc for PFHxS is often cited around 3.0. So, K_oc = 10^3 = 1000 L/kg.
k_oc_L_per_kg = 1000
print(f"\nAssuming K_oc for PFHxS = {k_oc_L_per_kg} L/kg (standard literature value).")

# Soil-water partition coefficient (Kd)
kd_L_per_kg = k_oc_L_per_kg * f_oc
print(f"Soil-Water Partition Coefficient (Kd) = {k_oc_L_per_kg} L/kg * {f_oc} = {kd_L_per_kg:.2f} L/kg")

# Convert bulk density to kg/L for unit consistency
# 1 m^3 = 1000 L
rho_b_kg_per_L = rho_b_kg_per_m3 / 1000
print(f"Bulk Density = {rho_b_kg_per_m3} kg/m^3 = {rho_b_kg_per_L} kg/L")

# Concentration in soil solution (C_sw)
# C_sw = C_soil / (Kd + (theta_w / rho_b))
c_sw_denominator = kd_L_per_kg + (theta_w / rho_b_kg_per_L)
c_sw_ug_per_L = c_soil_ug_per_kg / c_sw_denominator
print(f"Soil Solution Concentration (C_sw) = {c_soil_ug_per_kg:.4f} ug/kg / ({kd_L_per_kg:.2f} L/kg + ({theta_w} / {rho_b_kg_per_L} kg/L)) = {c_sw_ug_per_L:.4f} ug/L")

# --- Step 3: Calculate PFHxS Concentration in Produce ---
print("\n--- Step 3: Calculate PFHxS Concentration in Produce ---")
# We use the TSCF model. C_plant = TSCF * C_sw.
# We assume 1L of plant tissue is approximately 1kg.
c_fruit_ug_per_kg = tscf_fruit * c_sw_ug_per_L
print(f"Concentration in Fruits = {tscf_fruit} * {c_sw_ug_per_L:.4f} ug/L = {c_fruit_ug_per_kg:.4f} ug/kg")
c_legume_ug_per_kg = tscf_legume * c_sw_ug_per_L
print(f"Concentration in Legumes = {tscf_legume} * {c_sw_ug_per_L:.4f} ug/L = {c_legume_ug_per_kg:.4f} ug/kg")

# --- Step 4: Calculate Daily Intake of PFHxS ---
print("\n--- Step 4: Calculate Daily Intake ---")
# Convert intake rates from g to kg
intake_fruit_kg_day = intake_fruit_g_day / 1000
intake_legume_kg_day = intake_legume_g_day / 1000

# Daily intake from fruits
intake_from_fruit_ug_day = c_fruit_ug_per_kg * intake_fruit_kg_day * baf_fruit
print(f"Intake from Fruits = {c_fruit_ug_per_kg:.4f} ug/kg * {intake_fruit_kg_day} kg/day * {baf_fruit} = {intake_from_fruit_ug_day:.4f} ug/day")

# Daily intake from legumes
intake_from_legume_ug_day = c_legume_ug_per_kg * intake_legume_kg_day * baf_legume
print(f"Intake from Legumes = {c_legume_ug_per_kg:.4f} ug/kg * {intake_legume_kg_day} kg/day * {baf_legume} = {intake_from_legume_ug_day:.4f} ug/day")

# Total daily intake (TDI)
tdi_ug_day = intake_from_fruit_ug_day + intake_from_legume_ug_day
print(f"Total Daily Intake (TDI) = {intake_from_fruit_ug_day:.4f} ug/day + {intake_from_legume_ug_day:.4f} ug/day = {tdi_ug_day:.4f} ug/day")

# --- Step 5: Calculate Hazard Quotient (HQ) ---
print("\n--- Step 5: Calculate Hazard Quotient (HQ) ---")
# Chronic Daily Intake (CDI) normalized by body weight
cdi = tdi_ug_day / body_weight_kg
print(f"Chronic Daily Intake (CDI) = {tdi_ug_day:.4f} ug/day / {body_weight_kg} kg = {cdi:.6f} ug/kg/day")

# Hazard Quotient (HQ)
hq = cdi / ref_dose_ug_per_kg_day
print(f"\nHazard Quotient (HQ) = CDI / Reference Dose")
print(f"HQ = {cdi:.6f} ug/kg/day / {ref_dose_ug_per_kg_day} ug/kg/day = {hq:.4f}")

# Final Answer
final_answer = round(hq, 4)
