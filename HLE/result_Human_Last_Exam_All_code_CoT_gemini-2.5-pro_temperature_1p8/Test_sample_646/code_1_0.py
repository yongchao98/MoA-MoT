import math

# --- Given Data ---

# Contaminant Information
c_foam_ug_per_l = 1000000  # PFHxS concentration in foam (μg/L)
v_foam_l = 1000            # Volume of firefighting foam (L)

# Site Characteristics
area_m2 = 250000           # Area of land (m²)
depth_m = 0.6              # Depth of contamination (m)
bulk_density_kg_per_m3 = 1500 # Soil bulk density (kg/m³)

# Exposure Information
body_weight_kg = 80        # Body weight (kg)
ref_dose_ug_per_kg_day = 0.02 # Reference Dose for PFHxS (μg/kg/day)

# Consumption Data - Fruits
intake_fruit_g_day = 300
bioavailability_fruit = 0.5
puf_fruit = 0.1 # Plant Uptake Factor for fruits

# Consumption Data - Legumes
intake_legume_g_day = 50
bioavailability_legume = 0.3
puf_legume = 0.2 # Plant Uptake Factor for legumes

# --- Step-by-Step Calculation ---

# 1. Calculate the total mass of PFHxS applied
total_mass_pfhxs_ug = c_foam_ug_per_l * v_foam_l
print(f"Step 1: Calculate total mass of PFHxS applied")
print(f"  {c_foam_ug_per_l:,.0f} μg/L * {v_foam_l} L = {total_mass_pfhxs_ug:,.0f} μg\n")

# 2. Calculate the total mass of contaminated soil
volume_soil_m3 = area_m2 * depth_m
total_mass_soil_kg = volume_soil_m3 * bulk_density_kg_per_m3
print(f"Step 2: Calculate total mass of contaminated soil")
print(f"  Soil Volume = {area_m2:,.0f} m² * {depth_m} m = {volume_soil_m3:,.0f} m³")
print(f"  Soil Mass = {volume_soil_m3:,.0f} m³ * {bulk_density_kg_per_m3} kg/m³ = {total_mass_soil_kg:,.0f} kg\n")

# 3. Calculate the concentration of PFHxS in the soil (C_soil)
c_soil_ug_per_kg = total_mass_pfhxs_ug / total_mass_soil_kg
print(f"Step 3: Calculate PFHxS concentration in soil")
print(f"  C_soil = {total_mass_pfhxs_ug:,.0f} μg / {total_mass_soil_kg:,.0f} kg = {c_soil_ug_per_kg:.4f} μg/kg\n")

# 4. Calculate the concentration of PFHxS in the consumed produce
c_fruit_ug_per_kg = c_soil_ug_per_kg * puf_fruit
c_legume_ug_per_kg = c_soil_ug_per_kg * puf_legume
print(f"Step 4: Calculate PFHxS concentration in produce")
print(f"  C_fruit = {c_soil_ug_per_kg:.4f} μg/kg * {puf_fruit} = {c_fruit_ug_per_kg:.4f} μg/kg")
print(f"  C_legume = {c_soil_ug_per_kg:.4f} μg/kg * {puf_legume} = {c_legume_ug_per_kg:.4f} μg/kg\n")

# 5. Calculate the Total Daily Intake (TDI) of PFHxS
# Convert intake from g/day to kg/day
intake_fruit_kg_day = intake_fruit_g_day / 1000
intake_legume_kg_day = intake_legume_g_day / 1000

# Calculate intake from each source (µg/day)
daily_intake_fruit_ug = c_fruit_ug_per_kg * intake_fruit_kg_day * bioavailability_fruit
daily_intake_legume_ug = c_legume_ug_per_kg * intake_legume_kg_day * bioavailability_legume

# Calculate total daily intake and normalize by body weight
total_daily_intake_ug = daily_intake_fruit_ug + daily_intake_legume_ug
tdi_ug_per_kg_day = total_daily_intake_ug / body_weight_kg

print(f"Step 5: Calculate the Total Daily Intake (TDI)")
print(f"  Intake from fruits = {c_fruit_ug_per_kg:.4f} μg/kg * {intake_fruit_kg_day} kg/day * {bioavailability_fruit} = {daily_intake_fruit_ug:.4f} μg/day")
print(f"  Intake from legumes = {c_legume_ug_per_kg:.4f} μg/kg * {intake_legume_kg_day} kg/day * {bioavailability_legume} = {daily_intake_legume_ug:.4f} μg/day")
print(f"  Total intake = {daily_intake_fruit_ug:.4f} μg/day + {daily_intake_legume_ug:.4f} μg/day = {total_daily_intake_ug:.4f} μg/day")
print(f"  TDI = {total_daily_intake_ug:.4f} μg/day / {body_weight_kg} kg = {tdi_ug_per_kg_day:.6f} μg/kg/day\n")

# 6. Calculate the Hazard Quotient (HQ)
hazard_quotient = tdi_ug_per_kg_day / ref_dose_ug_per_kg_day
print(f"Step 6: Calculate the final Hazard Quotient (HQ)")
print(f"  The equation for the Hazard Quotient is:")
print(f"  HQ = Total Daily Intake (TDI) / Reference Dose (RfD)")
print(f"  HQ = {tdi_ug_per_kg_day:.6f} μg/kg/day / {ref_dose_ug_per_kg_day} μg/kg/day")
print(f"\nThe calculated Hazard Quotient is: {hazard_quotient:.2f}")
