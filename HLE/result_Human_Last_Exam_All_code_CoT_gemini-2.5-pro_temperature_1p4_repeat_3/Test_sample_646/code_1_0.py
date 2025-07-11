import math

# --- Given Data ---
# Soil properties
area = 250000  # m^2
depth = 0.6  # m
organic_carbon_content = 0.03  # 3% as a fraction
volumetric_water_content = 0.35  # L water/L soil
bulk_density_kg_m3 = 1500  # kg/m^3

# Contamination details
foam_volume = 1000  # L
pfhxs_conc_in_foam = 1000000  # ug/L

# Exposure details (Male)
body_weight = 80  # kg
reference_dose = 0.02  # ug/kg body weight per day

# Fruit consumption
fruit_intake_g = 300  # g/day
fruit_bioavailability = 0.5
fruit_puf = 0.1  # plant uptake factor
fruit_tscf = 5  # transpiration stream concentration factor

# Legume consumption
legume_intake_g = 50  # g/day
legume_bioavailability = 0.3
legume_puf = 0.2  # plant uptake factor
legume_tscf = 5  # transpiration stream concentration factor

# --- Calculations ---

# Step 1: Calculate the concentration of PFHxS in the soil (Cs)
# Total mass of PFHxS applied
total_pfhxs_mass_ug = foam_volume * pfhxs_conc_in_foam

# Total volume of affected soil
total_soil_volume_m3 = area * depth

# Total mass of affected soil
total_soil_mass_kg = total_soil_volume_m3 * bulk_density_kg_m3

# Concentration in soil (Cs) in ug/kg
cs_ug_per_kg = total_pfhxs_mass_ug / total_soil_mass_kg

# Step 2: Calculate the concentration of PFHxS in the soil solution (Cw)
# Convert bulk density to kg/L for unit consistency (1 m^3 = 1000 L)
bulk_density_kg_L = bulk_density_kg_m3 / 1000

# Assumption: The organic carbon-water partition coefficient (Koc) for PFHxS is not provided.
# A scientifically-based literature value is required. We will use a typical value where log(Koc) = 3.2.
log_koc = 3.2
koc_L_per_kg = 10**log_koc

# Soil partition coefficient (Kd) = Koc * fraction_of_organic_carbon
kd_L_per_kg = koc_L_per_kg * organic_carbon_content

# Concentration in soil water (Cw) using the equilibrium partitioning model: Cw = Cs / (Kd + (theta_w / rho_b))
cw_ug_per_L = cs_ug_per_kg / (kd_L_per_kg + (volumetric_water_content / bulk_density_kg_L))

# Step 3: Calculate PFHxS concentration in produce
# Assuming the density of plant tissue is ~1 kg/L, so ug/L is equivalent to ug/kg fresh weight.
c_fruit_ug_per_kg = cw_ug_per_L * fruit_tscf * fruit_puf
c_legume_ug_per_kg = cw_ug_per_L * legume_tscf * legume_puf

# Step 4: Calculate Chronic Daily Intake (CDI)
# Convert daily intake from g/day to kg/day
fruit_intake_kg = fruit_intake_g / 1000
legume_intake_kg = legume_intake_g / 1000

# CDI (ug/kg/day) = (Concentration * IntakeRate * Bioavailability) / BodyWeight
cdi_fruit = (c_fruit_ug_per_kg * fruit_intake_kg * fruit_bioavailability) / body_weight
cdi_legume = (c_legume_ug_per_kg * legume_intake_kg * legume_bioavailability) / body_weight
cdi_total = cdi_fruit + cdi_legume

# Step 5: Calculate the Hazard Quotient (HQ)
hq = cdi_total / reference_dose

# --- Output Results ---
print("This calculation determines the Hazard Quotient (HQ) for PFHxS exposure.")
print("NOTE: An assumption for the organic carbon-water partition coefficient (Koc) is required.")
print("A log(Koc) of 3.2, a standard literature value for PFHxS, has been used.")
print("\n--- Intermediate Calculations ---")
print(f"1. Concentration in Soil (Cs): {cs_ug_per_kg:.4f} µg/kg")
print(f"2. Concentration in Soil Water (Cw): {cw_ug_per_L:.4f} µg/L")
print(f"3. Concentration in Fruits (C_fruit): {c_fruit_ug_per_kg:.4f} µg/kg")
print(f"4. Concentration in Legumes (C_legume): {c_legume_ug_per_kg:.4f} µg/kg")
print(f"5. Total Daily Intake (CDI_total): {cdi_total:.8f} µg/kg/day")

print("\n--- Final Hazard Quotient Calculation ---")
print("Hazard Quotient (HQ) = Total Daily Intake / Reference Dose")
print("Equation: HQ = [ (C_fruit * IR_fruit * BF_fruit / BW) + (C_legume * IR_legume * BF_legume / BW) ] / RfD")
print(f"HQ = [ ({c_fruit_ug_per_kg:.4f} * {fruit_intake_kg} * {fruit_bioavailability} / {body_weight}) + ({c_legume_ug_per_kg:.4f} * {legume_intake_kg} * {legume_bioavailability} / {body_weight}) ] / {reference_dose}")
print(f"HQ = ({cdi_total:.8f}) / ({reference_dose})")
print(f"Final Hazard Quotient = {hq:.6f}")
