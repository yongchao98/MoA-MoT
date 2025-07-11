import math

# --- Given Parameters ---

# Soil and Contamination
area_m2 = 250000
depth_m = 0.6
foc = 0.03  # Fraction of organic carbon (3%)
theta_w_ratio = 0.35  # Volumetric water content (L water/L soil)
bulk_density_kg_m3 = 1500
foam_volume_L = 1000
pfhxs_conc_in_foam_ug_L = 1000000

# Exposure Parameters
body_weight_kg = 80
# Fruits
intake_rate_fruit_g_day = 300
tscf_fruit = 5
# Legumes
intake_rate_legume_g_day = 50
tscf_legume = 5

# Toxicological Data
ref_dose_ug_kg_day = 0.02

# --- Assumptions ---
# Koc for PFHxS is not given, a literature value of 1000 L/kg (log Koc = 3.0) is assumed.
koc_L_kg = 1000
# The 'plant uptake factor' and 'bioavailability factor' are not used in this model.
# The Hazard Quotient calculation compares total intake to the RfD, not the absorbed dose.
# The concentration in the edible part of the plant is assumed to be equal to the concentration
# in the transpiration stream (plant sap), and the density of the food is ~1 kg/L.

# --- Step 1: Calculate Total PFHxS Mass ---
total_pfhxs_ug = foam_volume_L * pfhxs_conc_in_foam_ug_L
print(f"Step 1: Total PFHxS applied = {foam_volume_L} L * {pfhxs_conc_in_foam_ug_L} µg/L = {total_pfhxs_ug:,.0f} µg")

# --- Step 2: Calculate Total Soil Mass ---
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * bulk_density_kg_m3
print(f"Step 2: Total soil mass = ({area_m2} m² * {depth_m} m) * {bulk_density_kg_m3} kg/m³ = {soil_mass_kg:,.0f} kg")

# --- Step 3: Calculate PFHxS Concentration in Soil (Cs) ---
cs_ug_kg = total_pfhxs_ug / soil_mass_kg
print(f"Step 3: Concentration in soil (Cs) = {total_pfhxs_ug:,.0f} µg / {soil_mass_kg:,.0f} kg = {cs_ug_kg:.4f} µg/kg")

# --- Step 4: Calculate PFHxS Concentration in Soil Solution (Cw) ---
# Calculate soil-water partition coefficient (Kd)
kd_L_kg = koc_L_kg * foc
print(f"Step 4a: Soil-water partition coefficient (Kd) = {koc_L_kg} L/kg * {foc} = {kd_L_kg:.1f} L/kg")

# Convert theta_w from ratio to absolute volume (L/m³)
theta_w_abs_L_m3 = theta_w_ratio * 1000  # 1 m³ = 1000 L

# Calculate Cw using the partitioning formula: Cw = (Cs * ρb) / (θw_abs + Kd * ρb)
# Units: Cw[µg/L] = (Cs[µg/kg] * ρb[kg/m³]) / (θw[L/m³] + Kd[L/kg] * ρb[kg/m³])
cw_ug_L = (cs_ug_kg * bulk_density_kg_m3) / (theta_w_abs_L_m3 + kd_L_kg * bulk_density_kg_m3)
print(f"Step 4b: Concentration in soil solution (Cw) = {cw_ug_L:.4f} µg/L")

# --- Step 5: Calculate PFHxS Concentration in Food (C_food) ---
# Assuming C_food = Cw * TSCF and density of food is ~1 kg/L
c_food_ug_L = cw_ug_L * tscf_fruit # TSCF is the same for both
c_food_ug_kg = c_food_ug_L
c_food_ug_g = c_food_ug_kg / 1000
print(f"Step 5: Concentration in food (C_food) = {cw_ug_L:.4f} µg/L * {tscf_fruit} = {c_food_ug_L:.4f} µg/L ≈ {c_food_ug_g:.8f} µg/g")

# --- Step 6: Calculate Total Daily Intake ---
daily_intake_fruit_ug = c_food_ug_g * intake_rate_fruit_g_day
daily_intake_legume_ug = c_food_ug_g * intake_rate_legume_g_day
total_daily_intake_ug_day = daily_intake_fruit_ug + daily_intake_legume_ug
print(f"Step 6: Total daily intake = ({c_food_ug_g:.8f} µg/g * {intake_rate_fruit_g_day} g/day) + ({c_food_ug_g:.8f} µg/g * {intake_rate_legume_g_day} g/day) = {total_daily_intake_ug_day:.4f} µg/day")

# --- Step 7: Calculate Chronic Daily Intake (CDI) ---
cdi_ug_kg_day = total_daily_intake_ug_day / body_weight_kg
print(f"Step 7: Chronic Daily Intake (CDI) = {total_daily_intake_ug_day:.4f} µg/day / {body_weight_kg} kg = {cdi_ug_kg_day:.8f} µg/kg/day")

# --- Step 8: Calculate Hazard Quotient (HQ) ---
hazard_quotient = cdi_ug_kg_day / ref_dose_ug_kg_day
print(f"Step 8: Hazard Quotient (HQ) = {cdi_ug_kg_day:.8f} µg/kg/day / {ref_dose_ug_kg_day} µg/kg/day = {hazard_quotient:.4f}")

print("\n--- Final Hazard Quotient Calculation ---")
final_equation = (f"HQ = [ (C_food * IR_fruit) + (C_food * IR_legume) ] / BW / RfD\n"
                  f"HQ = [ ({c_food_ug_g:.8f} µg/g * {intake_rate_fruit_g_day} g/day) + ({c_food_ug_g:.8f} µg/g * {intake_rate_legume_g_day} g/day) ] / {body_weight_kg} kg / {ref_dose_ug_kg_day} µg/kg/day\n"
                  f"HQ = [ {daily_intake_fruit_ug:.4f} µg/day + {daily_intake_legume_ug:.4f} µg/day ] / {body_weight_kg} kg / {ref_dose_ug_kg_day} µg/kg/day\n"
                  f"HQ = {total_daily_intake_ug_day:.4f} µg/day / {body_weight_kg} kg / {ref_dose_ug_kg_day} µg/kg/day\n"
                  f"HQ = {cdi_ug_kg_day:.8f} µg/kg/day / {ref_dose_ug_kg_day} µg/kg/day\n"
                  f"HQ = {hazard_quotient:.4f}")
print(final_equation)

# Final answer in the required format
print(f"\n<<<{hazard_quotient:.4f}>>>")