import math

# --- Given Data ---

# Soil and Contamination Properties
area_m2 = 250000  # m^2
depth_m = 0.6  # m
f_oc = 0.03  # organic carbon content fraction
theta_w = 0.35  # volumetric water content (L water / L soil)
rho_b_kg_m3 = 1500  # bulk density (kg/m^3)
foam_volume_L = 1000  # L
foam_conc_ug_L = 1000000  # ug/L

# PFHxS Property (Standard literature value, as not provided in the problem)
# This is a key assumption.
K_oc_L_kg = 251  # Organic carbon-water partitioning coefficient (L/kg)

# Exposure Scenario
body_weight_kg = 80  # kg

# Fruit Consumption
IR_fruit_g_day = 300  # g/day
BF_fruit = 0.5  # bioavailability factor
TSCF_fruit = 5  # transpiration stream concentration factor

# Legume Consumption
IR_legume_g_day = 50  # g/day
BF_legume = 0.3  # bioavailability factor
TSCF_legume = 5  # transpiration stream concentration factor

# Toxicological Data
RfD_ug_kg_day = 0.02  # reference dose (ug/kg/day)

# --- Step-by-Step Calculation ---

# 1. Calculate total mass of PFHxS applied in micrograms (ug)
total_PFHxS_ug = foam_volume_L * foam_conc_ug_L

# 2. Calculate total volume and mass of contaminated soil
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * rho_b_kg_m3

# 3. Calculate concentration of PFHxS in bulk soil (Cs) in ug/kg
Cs_ug_kg = total_PFHxS_ug / soil_mass_kg

# 4. Calculate concentration in soil solution (Cw) in ug/L
# First, calculate soil-water partitioning coefficient (Kd)
Kd_L_kg = K_oc_L_kg * f_oc
# Then, calculate Cw using the partitioning equation: Cw = (Cs * rho_b) / (theta_w_converted + Kd * rho_b)
# Convert volumetric water content theta_w from L/L to L/m^3
theta_w_L_m3 = theta_w * 1000  # 1 m^3 = 1000 L
# Now calculate Cw
numerator_cw = Cs_ug_kg * rho_b_kg_m3
denominator_cw = theta_w_L_m3 + (Kd_L_kg * rho_b_kg_m3)
Cw_ug_L = numerator_cw / denominator_cw

# 5. Calculate PFHxS concentration in plants (C_plant) in ug/kg fresh weight
# We assume 1 L of plant matter has a mass of 1 kg.
C_plant_fruit_ug_kg = Cw_ug_L * TSCF_fruit
C_plant_legume_ug_kg = Cw_ug_L * TSCF_legume

# 6. Calculate Average Daily Dose (ADD) for each food in ug/kg/day
# Convert intake rates from g/day to kg/day
IR_fruit_kg_day = IR_fruit_g_day / 1000
IR_legume_kg_day = IR_legume_g_day / 1000

ADD_fruit_ug_kg_day = (C_plant_fruit_ug_kg * IR_fruit_kg_day * BF_fruit) / body_weight_kg
ADD_legume_ug_kg_day = (C_plant_legume_ug_kg * IR_legume_kg_day * BF_legume) / body_weight_kg

# Calculate total daily intake
Total_ADD_ug_kg_day = ADD_fruit_ug_kg_day + ADD_legume_ug_kg_day

# 7. Calculate Hazard Quotient (HQ)
HQ = Total_ADD_ug_kg_day / RfD_ug_kg_day

# --- Print Results ---
print("--- Calculation Steps ---")
print(f"1. Total PFHxS applied: {total_PFHxS_ug:,.0f} µg")
print(f"2. Total mass of contaminated soil: {soil_mass_kg:,.0f} kg")
print(f"3. Concentration in bulk soil (Cs): {Cs_ug_kg:.4f} µg/kg")
print(f"4. Concentration in soil solution (Cw): {Cw_ug_L:.4f} µg/L")
print(f"5. Concentration in fruits: {C_plant_fruit_ug_kg:.4f} µg/kg")
print(f"   Concentration in legumes: {C_plant_legume_ug_kg:.4f} µg/kg")
print(f"6. Total daily intake (ADD): {Total_ADD_ug_kg_day:.6f} µg/kg/day")
print("\n--- Final Hazard Quotient Calculation ---")
print("Hazard Quotient (HQ) = Total Daily Intake / Reference Dose")
print(f"HQ = {Total_ADD_ug_kg_day:.6f} µg/kg/day / {RfD_ug_kg_day} µg/kg/day")
print(f"HQ = {HQ:.3f}")

if HQ < 1:
    print("\nResult: The Hazard Quotient is less than 1, suggesting that adverse health effects are unlikely based on this exposure scenario.")
else:
    print("\nResult: The Hazard Quotient is 1 or greater, suggesting a potential health risk.")

# Final answer in the required format
# print(f"\n<<<{HQ}>>>")