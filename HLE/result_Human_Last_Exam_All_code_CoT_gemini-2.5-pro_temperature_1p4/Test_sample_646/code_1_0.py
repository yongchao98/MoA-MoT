import math

# --- Step 1: Define Given Data and Constants ---
# PFHxS and Foam Data
foam_volume_L = 1000  # L
pfhxs_conc_in_foam_ug_per_L = 1000000  # μg/L

# Soil Data
area_m2 = 250000  # m²
depth_m = 0.6  # m
f_oc = 0.03  # organic carbon content (3%)
theta_w = 0.35  # volumetric water content (L water/L soil)
rho_b_kg_per_m3 = 1500  # bulk density in kg/m³

# Human and Consumption Data
body_weight_kg = 80  # kg

# Fruits
intake_fruit_g_per_day = 300
bf_fruit = 0.5
puf_fruit = 0.1

# Legumes
intake_legume_g_per_day = 50
bf_legume = 0.3
puf_legume = 0.2

# General Factors
tscf = 5  # Transpiration stream concentration factor

# Toxicological Data
rfd_ug_per_kg_day = 0.02  # Reference Dose in μg/kg/day

# Assumption
# The organic carbon-water partition coefficient (Koc) is required.
# A literature value for log(Koc) for PFHxS is used.
log_koc_pfhxs = 3.15

# --- Step 2: Calculate Initial Contamination ---
print("--- Calculating Soil Contamination ---")
total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L
print(f"1. Total mass of PFHxS applied: {total_pfhxs_mass_ug:,.0f} μg")

soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * rho_b_kg_per_m3
print(f"2. Total mass of contaminated soil: {soil_mass_kg:,.0f} kg")

c_soil_ug_per_kg = total_pfhxs_mass_ug / soil_mass_kg
print(f"3. Concentration of PFHxS in bulk soil (C_soil): {c_soil_ug_per_kg:.4f} μg/kg")

# --- Step 3: Calculate Concentration in Soil Solution (Pore Water) ---
# Convert bulk density to kg/L for consistency in units (1 m³ = 1000 L)
rho_b_kg_per_L = rho_b_kg_per_m3 / 1000

# Calculate Koc and Kd
k_oc = 10**log_koc_pfhxs
kd = k_oc * f_oc

# Calculate concentration in pore water (C_w) using the partition equation:
# C_w = C_soil / (Kd + theta_w / rho_b)
c_w_ug_per_L = c_soil_ug_per_kg / (kd + (theta_w / rho_b_kg_per_L))
print(f"4. Concentration of PFHxS in soil solution (C_w): {c_w_ug_per_L:.4f} μg/L\n")

# --- Step 4: Calculate Concentration in Produce ---
# Assuming density of produce is ~1 kg/L, so C_plant in μg/kg is numerically equal to C_plant in μg/L.
print("--- Calculating Contaminant Levels in Food ---")
c_fruit_ug_per_kg = c_w_ug_per_L * tscf * puf_fruit
print(f"5. Concentration in fruits: {c_fruit_ug_per_kg:.4f} μg/kg")

c_legume_ug_per_kg = c_w_ug_per_L * tscf * puf_legume
print(f"6. Concentration in legumes: {c_legume_ug_per_kg:.4f} μg/kg\n")

# --- Step 5: Calculate Daily Intake ---
print("--- Calculating Daily Intake ---")
# Convert intake rates from g/day to kg/day
intake_fruit_kg_per_day = intake_fruit_g_per_day / 1000
intake_legume_kg_per_day = intake_legume_g_per_day / 1000

# Calculate intake for each food type, normalized by body weight
add_fruit = (c_fruit_ug_per_kg * intake_fruit_kg_per_day * bf_fruit) / body_weight_kg
print(f"7. Daily intake from fruits (normalized): {add_fruit:.8f} μg/kg/day")

add_legume = (c_legume_ug_per_kg * intake_legume_kg_per_day * bf_legume) / body_weight_kg
print(f"8. Daily intake from legumes (normalized): {add_legume:.8f} μg/kg/day")

# Calculate Total Daily Intake (Total ADD)
total_add = add_fruit + add_legume
print(f"9. Total Daily Intake (ADD): {total_add:.8f} μg/kg/day\n")

# --- Step 6: Calculate Hazard Quotient ---
print("--- Calculating Hazard Quotient (HQ) ---")
# HQ = Total Daily Intake / Reference Dose
hq = total_add / rfd_ug_per_kg_day

print("The Hazard Quotient is calculated as:")
print(f"HQ = Total Daily Intake (ADD) / Reference Dose (RfD)")
print(f"HQ = {total_add:.8f} μg/kg/day / {rfd_ug_per_kg_day} μg/kg/day")
print(f"Final Hazard Quotient (HQ): {hq:.4f}")

<<<0.0059>>>