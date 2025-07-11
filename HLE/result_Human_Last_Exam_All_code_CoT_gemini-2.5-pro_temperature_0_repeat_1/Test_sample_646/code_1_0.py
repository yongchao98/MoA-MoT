import math

# Step 1: Define constants and given variables

# Contamination Data
foam_volume_L = 1000
pfhxs_conc_in_foam_ug_per_L = 1000000

# Soil Data
area_m2 = 250000
depth_m = 0.6
f_oc = 0.03  # Organic carbon content (3%)
theta_w_L_per_L = 0.35  # Volumetric water content
rho_b_kg_per_m3 = 1500  # Bulk density

# Exposure Data (Man)
body_weight_kg = 80

# Fruit Consumption Data
fruit_intake_g_per_day = 300
fruit_bf = 0.5  # Bioavailability factor
fruit_puf = 0.1  # Plant uptake factor
fruit_tscf = 5  # Transpiration stream concentration factor

# Legume Consumption Data
legume_intake_g_per_day = 50
legume_bf = 0.3  # Bioavailability factor
legume_puf = 0.2  # Plant uptake factor
legume_tscf = 5  # Transpiration stream concentration factor

# Toxicological Data
rfd_ug_per_kg_day = 0.02  # Reference Dose

# Assumed value for PFHxS partitioning
# K_oc (organic carbon-water partition coefficient) is not given.
# A common literature value for log(K_oc) is 2.8.
log_koc = 2.8
k_oc_L_per_kg = 10**log_koc

# Step 2: Calculate total PFHxS mass and soil mass
total_pfhxs_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L
soil_volume_m3 = area_m2 * depth_m
total_soil_mass_kg = soil_volume_m3 * rho_b_kg_per_m3

# Step 3: Calculate concentration in soil (C_soil)
c_soil_ug_per_kg = total_pfhxs_ug / total_soil_mass_kg

# Step 4: Calculate concentration in soil solution (C_solution)
kd_L_per_kg = k_oc_L_per_kg * f_oc
# Convert volumetric water content from L/L to L/m^3
theta_w_L_per_m3 = theta_w_L_per_L * 1000
# Formula: C_solution = (C_soil * rho_b) / (Kd * rho_b + theta_w)
c_solution_ug_per_L = (c_soil_ug_per_kg * rho_b_kg_per_m3) / (kd_L_per_kg * rho_b_kg_per_m3 + theta_w_L_per_m3)

# Step 5: Calculate total daily intake (TDI)
# Convert daily intake from g/day to kg/day
fruit_intake_kg_per_day = fruit_intake_g_per_day / 1000
legume_intake_kg_per_day = legume_intake_g_per_day / 1000

# Calculate concentration in fruits and legumes (assuming plant tissue density is 1 kg/L)
c_fruit_ug_per_kg = c_solution_ug_per_L * fruit_tscf * fruit_puf
c_legume_ug_per_kg = c_solution_ug_per_L * legume_tscf * legume_puf

# Calculate absorbed daily dose from each food source
# Formula: (Concentration_in_food * Intake_Rate * Bioavailability_Factor) / Body_Weight
di_fruit_ug_per_kg_day = (c_fruit_ug_per_kg * fruit_intake_kg_per_day * fruit_bf) / body_weight_kg
di_legume_ug_per_kg_day = (c_legume_ug_per_kg * legume_intake_kg_per_day * legume_bf) / body_weight_kg

# Calculate Total Daily Intake (TDI)
tdi_ug_per_kg_day = di_fruit_ug_per_kg_day + di_legume_ug_per_kg_day

# Step 6: Calculate the Hazard Quotient (HQ)
# Formula: HQ = TDI / RfD
hazard_quotient = tdi_ug_per_kg_day / rfd_ug_per_kg_day

# Step 7: Print the final calculation and result
print("Calculation of the Hazard Quotient (HQ):")
print("-" * 40)
print(f"Total Daily Intake (TDI) = {tdi_ug_per_kg_day:.9f} µg/kg/day")
print(f"Reference Dose (RfD) = {rfd_ug_per_kg_day} µg/kg/day")
print("\nHazard Quotient (HQ) = Total Daily Intake / Reference Dose")
print(f"HQ = {tdi_ug_per_kg_day:.9f} / {rfd_ug_per_kg_day}")
print(f"HQ = {hazard_quotient:.5f}")

<<<0.01305>>>