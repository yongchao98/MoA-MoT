import math

# --- Given Data ---

# Soil and Contamination Properties
area_m2 = 250000
depth_m = 0.6
f_oc = 0.03  # organic carbon content (3%)
theta_w = 0.35  # volumetric water content (L/L)
rho_b_m3 = 1500  # bulk density (kg/m^3)

# Foam and PFHxS Properties
foam_volume_L = 1000
pfhxs_conc_in_foam_ug_per_L = 1000000

# Man and Consumption Data
body_weight_kg = 80
rfd_ug_per_kg_day = 0.02  # Reference Dose

# Fruit Data
ir_fruit_g = 300
baf_fruit = 0.5
puf_fruit = 0.1
tscf = 5  # Transpiration stream concentration factor (same for both)

# Legume Data
ir_legume_g = 50
baf_legume = 0.3
puf_legume = 0.2

# --- Calculations ---

# Step 1: Calculate initial soil concentration (C_soil)
total_pfhxs_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * rho_b_m3
c_soil_ug_per_kg = total_pfhxs_ug / soil_mass_kg

# Step 2: Calculate concentration in soil water (Cw)
# Assume Koc for PFHxS (not provided in the problem)
k_oc_L_per_kg = 1585  # Literature value for PFHxS (log Koc = 3.2)
k_d_L_per_kg = k_oc_L_per_kg * f_oc

# Convert bulk density to kg/L for unit consistency in the formula
rho_b_L = rho_b_m3 / 1000.0  # 1 m^3 = 1000 L

# Cw = C_soil / (Kd + (theta_w / rho_b))
c_w_ug_per_L = c_soil_ug_per_kg / (k_d_L_per_kg + (theta_w / rho_b_L))

# Step 3: Calculate concentration in plants (fruits and legumes)
# Assuming 1L of plant sap is equivalent to 1kg of plant fresh weight
c_fruit_ug_per_kg = c_w_ug_per_L * tscf * puf_fruit
c_legume_ug_per_kg = c_w_ug_per_L * tscf * puf_legume

# Step 4: Calculate total daily intake (TDI)
# Convert intake rates from g/day to kg/day
ir_fruit_kg = ir_fruit_g / 1000.0
ir_legume_kg = ir_legume_g / 1000.0

# Absorbed dose from each food source (ug/day)
absorbed_dose_fruit_ug_day = c_fruit_ug_per_kg * ir_fruit_kg * baf_fruit
absorbed_dose_legume_ug_day = c_legume_ug_per_kg * ir_legume_kg * baf_legume

# Total daily intake (TDI) normalized by body weight (ug/kg/day)
tdi_ug_per_kg_day = (absorbed_dose_fruit_ug_day + absorbed_dose_legume_ug_day) / body_weight_kg

# Step 5: Calculate Hazard Quotient (HQ)
hq = tdi_ug_per_kg_day / rfd_ug_per_kg_day

# --- Output the Final Equation and Result ---
print("This calculation determines the Hazard Quotient (HQ) based on the man's daily consumption.")
print("\nFirst, the daily intake from each food source is calculated and combined:")
print(f"Intake from Fruits (µg/day) = Concentration ({c_fruit_ug_per_kg:.6f} µg/kg) * Intake Rate ({ir_fruit_kg} kg/day) * Bioavailability ({baf_fruit}) = {absorbed_dose_fruit_ug_day:.6f}")
print(f"Intake from Legumes (µg/day) = Concentration ({c_legume_ug_per_kg:.6f} µg/kg) * Intake Rate ({ir_legume_kg} kg/day) * Bioavailability ({baf_legume}) = {absorbed_dose_legume_ug_day:.6f}")
print("\nNext, the total daily intake is normalized by body weight and compared to the reference dose to find the HQ.")
print("\n--- Final Equation ---")
print(f"Hazard Quotient = (Intake from Fruits + Intake from Legumes) / Body Weight / Reference Dose")
print(f"Hazard Quotient = ({absorbed_dose_fruit_ug_day:.6f} µg/day + {absorbed_dose_legume_ug_day:.6f} µg/day) / {body_weight_kg} kg / {rfd_ug_per_kg_day} µg/kg/day")
print(f"Hazard Quotient = {hq:.6f}")
