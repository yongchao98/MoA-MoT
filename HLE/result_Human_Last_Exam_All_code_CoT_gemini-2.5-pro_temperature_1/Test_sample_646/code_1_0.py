import math

# Plan:
# 1. Define all given parameters.
# 2. Calculate the total concentration of PFHxS in soil (Ct).
# 3. Calculate the concentration of PFHxS in the soil water/solution (Cw).
# 4. Calculate the concentration of PFHxS in the consumed produce (fruits and legumes).
# 5. Calculate the Chronic Daily Intake (CDI) from both food sources.
# 6. Calculate the final Hazard Quotient (HQ).

# --- Step 1: Define Given Parameters ---

# Soil and Contaminant Properties
soil_area_m2 = 250000
soil_depth_m = 0.6
organic_carbon_fraction = 0.03  # 3%
volumetric_water_content = 0.35  # L water / L soil
bulk_density_kg_m3 = 1500
foam_volume_L = 1000
pfhxs_in_foam_ug_L = 1000000

# Assumption: Koc for PFHxS is not provided. A literature value will be used.
# logKoc for PFHxS is ~2.6, so Koc is 10^2.6.
koc_L_kg = 10**2.6 # ~398.1 L/kg

# Exposure Scenario
body_weight_kg = 80
rfd_ug_kg_day = 0.02  # Reference Dose

# Fruits Consumption
fruit_intake_g_day = 300
fruit_bioavailability = 0.5
fruit_tscf = 5  # Transpiration Stream Concentration Factor

# Legumes Consumption
legume_intake_g_day = 50
legume_bioavailability = 0.3
legume_tscf = 5  # Transpiration Stream Concentration Factor

print("--- Hazard Quotient Calculation for PFHxS Exposure ---")
print(f"Assumption: The organic carbon-water partitioning coefficient (Koc) for PFHxS is {koc_L_kg:.2f} L/kg.\n")

# --- Step 2: Calculate Total Concentration in Soil (Ct) ---

# Total mass of PFHxS added to the soil
total_pfhxs_ug = foam_volume_L * pfhxs_in_foam_ug_L

# Total volume of contaminated soil
soil_volume_m3 = soil_area_m2 * soil_depth_m

# Total mass of contaminated soil
total_soil_mass_kg = soil_volume_m3 * bulk_density_kg_m3

# Total concentration (Ct) in soil
ct_ug_kg = total_pfhxs_ug / total_soil_mass_kg
print(f"Step 1: Total PFHxS concentration in soil (Ct) = {ct_ug_kg:.4f} µg/kg")

# --- Step 3: Calculate Concentration in Soil Solution (Cw) ---

# Distribution coefficient (Kd)
kd_L_kg = koc_L_kg * organic_carbon_fraction

# Convert bulk density to kg/L for unit consistency (1 m^3 = 1000 L)
bulk_density_kg_L = bulk_density_kg_m3 / 1000

# Concentration in soil water (Cw)
# Formula: Cw = Ct / (Kd + θw / ρb)
cw_ug_L = ct_ug_kg / (kd_L_kg + (volumetric_water_content / bulk_density_kg_L))
print(f"Step 2: PFHxS concentration in soil solution (Cw) = {cw_ug_L:.4f} µg/L")

# --- Step 4: Calculate Concentration in Plants (C_plant) ---

# C_plant = Cw * TSCF. Assuming plant density is ~1 kg/L, so C_plant is in µg/kg.
c_fruit_ug_kg = cw_ug_L * fruit_tscf
c_legume_ug_kg = cw_ug_L * legume_tscf
print(f"Step 3: PFHxS concentration in fruits = {c_fruit_ug_kg:.4f} µg/kg")
print(f"         PFHxS concentration in legumes = {c_legume_ug_kg:.4f} µg/kg")


# --- Step 5: Calculate Chronic Daily Intake (CDI) ---

# Convert intake rates from g/day to kg/day
fruit_intake_kg_day = fruit_intake_g_day / 1000
legume_intake_kg_day = legume_intake_g_day / 1000

# CDI for fruits
cdi_fruit = (c_fruit_ug_kg * fruit_intake_kg_day * fruit_bioavailability) / body_weight_kg

# CDI for legumes
cdi_legume = (c_legume_ug_kg * legume_intake_kg_day * legume_bioavailability) / body_weight_kg

# Total CDI
cdi_total = cdi_fruit + cdi_legume
print(f"\nStep 4: Total Chronic Daily Intake (CDI) = {cdi_total:.6f} µg/kg/day")

# --- Step 6: Calculate Hazard Quotient (HQ) ---

# HQ = Total CDI / Reference Dose (RfD)
hq = cdi_total / rfd_ug_kg_day

print("\n--- Final Result ---")
print("The Hazard Quotient (HQ) is the ratio of the Chronic Daily Intake to the Reference Dose.")
# Final equation output as requested
print(f"Hazard Quotient (HQ) = (Total CDI) / (Reference Dose)")
print(f"Hazard Quotient (HQ) = {cdi_total:.6f} µg/kg/day / {rfd_ug_kg_day} µg/kg/day")
print(f"The calculated Hazard Quotient is: {hq:.4f}")
