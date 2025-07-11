import math

# --- Given Data ---

# Contamination Data
foam_volume_L = 1000
pfhxs_conc_in_foam_ug_per_L = 1000000

# Soil Properties
area_m2 = 250000
depth_m = 0.6
organic_carbon_content_pct = 0.03  # 3% as a fraction
volumetric_water_content = 0.35  # theta_w, L water/L soil
bulk_density_kg_per_m3 = 1500

# Exposure Scenario
body_weight_kg = 80

# Fruit Consumption
ir_fruit_g_day = 300
bf_fruit = 0.5
puf_fruit = 0.1
tscf_fruit = 5

# Legume Consumption
ir_legume_g_day = 50
bf_legume = 0.3
puf_legume = 0.2
tscf_legume = 5

# Toxicological Data
rfd_ug_per_kg_day = 0.02

# --- Calculations ---

# Step 1: Calculate total PFHxS mass applied
total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_per_L

# Step 2: Calculate soil volume and mass
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * bulk_density_kg_per_m3

# Step 3: Calculate PFHxS concentration in soil (Cs)
# This is the mass of PFHxS per kg of dry soil
cs_ug_per_kg = total_pfhxs_mass_ug / soil_mass_kg

# Step 4: Calculate PFHxS concentration in soil water (Cw)
# First, estimate the soil-water partition coefficient (Kd)
log_kow = 4.5  # Assumed log Octanol-Water Partition Coefficient for PFHxS
log_koc = 0.52 * log_kow + 1.02 # Empirical estimation of Koc from Kow
koc = 10**log_koc
kd_L_per_kg = koc * organic_carbon_content_pct

# Convert bulk density to kg/L
bulk_density_kg_per_L = bulk_density_kg_per_m3 / 1000

# Calculate Cw using equilibrium partitioning model
# Cw = (Total mass / Total volume) / (water_content + Kd * bulk_density)
total_soil_volume_L = soil_volume_m3 * 1000
total_conc_ug_per_L = total_pfhxs_mass_ug / total_soil_volume_L
cw_ug_per_L = total_conc_ug_per_L / (volumetric_water_content + kd_L_per_kg * bulk_density_kg_per_L)

# Step 5: Calculate PFHxS concentration in plants (Cplant)
c_plant_fruit_ug_per_kg = cw_ug_per_L * tscf_fruit * puf_fruit
c_plant_legume_ug_per_kg = cw_ug_per_L * tscf_legume * puf_legume

# Step 6: Calculate daily intake (DI) from each food source
# Convert intake rates from g/day to kg/day
ir_fruit_kg_day = ir_fruit_g_day / 1000
ir_legume_kg_day = ir_legume_g_day / 1000

di_fruit_ug_day = c_plant_fruit_ug_per_kg * ir_fruit_kg_day * bf_fruit
di_legume_ug_day = c_plant_legume_ug_per_kg * ir_legume_kg_day * bf_legume

# Step 7: Calculate Chronic Daily Intake (CDI)
total_daily_intake_ug_day = di_fruit_ug_day + di_legume_ug_day
cdi_ug_per_kg_day = total_daily_intake_ug_day / body_weight_kg

# Step 8: Calculate Hazard Quotient (HQ)
hazard_quotient = cdi_ug_per_kg_day / rfd_ug_per_kg_day

# --- Print Results ---
print("--- Calculation Steps ---")
print(f"1. Total PFHxS Applied: {total_pfhxs_mass_ug:,.0f} µg")
print(f"2. Concentration in Soil (Cs): {cs_ug_per_kg:.4f} µg/kg")
print(f"3. Concentration in Soil Water (Cw): {cw_ug_per_L:.4f} µg/L")
print(f"4. Concentration in Fruit: {c_plant_fruit_ug_per_kg:.4f} µg/kg")
print(f"5. Concentration in Legumes: {c_plant_legume_ug_per_kg:.4f} µg/kg")
print(f"6. Total Daily Intake (TDI): {total_daily_intake_ug_day:.4f} µg/day")
print(f"7. Chronic Daily Intake (CDI): {cdi_ug_per_kg_day:.6f} µg/kg/day")
print("\n--- Final Hazard Quotient Calculation ---")
print(f"HQ = CDI / RfD")
print(f"HQ = {cdi_ug_per_kg_day:.6f} µg/kg/day / {rfd_ug_per_kg_day} µg/kg/day")
print(f"Hazard Quotient (HQ) = {hazard_quotient:.4f}")

# The final answer in the required format
print(f"\n<<<{hazard_quotient:.4f}>>>")