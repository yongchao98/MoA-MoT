import math

# --- Step 1: Define initial data and constants ---
# Soil and contamination properties
area = 250000  # m^2
depth = 0.6  # m
bulk_density_m3 = 1500  # kg/m^3
f_oc = 0.03  # organic carbon fraction
theta_w = 0.35  # volumetric water content (L water / L soil)

# Foam and PFHxS properties
foam_volume = 1000  # Litres
pfhxs_in_foam = 1000000  # ug/L
log_koc_pfhxs = 2.5 # Assumed standard value for PFHxS log Koc

# Human exposure properties
body_weight = 80  # kg
ref_dose = 0.02  # ug/kg body weight per day

# Fruit properties
intake_fruit_g = 300  # g/day
bf_fruit = 0.5  # bioavailability factor
puf_fruit = 0.1  # plant uptake factor
tscf = 5  # transpiration stream concentration factor

# Legume properties
intake_legume_g = 50  # g/day
bf_legume = 0.3  # bioavailability factor
puf_legume = 0.2  # plant uptake factor

# --- Step 2: Calculate total PFHxS mass and soil mass ---
total_pfhxs_mass = foam_volume * pfhxs_in_foam
soil_volume_m3 = area * depth
total_soil_mass_kg = soil_volume_m3 * bulk_density_m3

# --- Step 3: Calculate PFHxS concentration in soil (C_soil) ---
c_soil = total_pfhxs_mass / total_soil_mass_kg

# --- Step 4: Calculate PFHxS concentration in soil solution (C_solution) ---
# Koc from logKoc
koc = 10**log_koc_pfhxs
# Kd from Koc
kd = koc * f_oc
# Convert bulk density to kg/L for unit consistency
bulk_density_L = bulk_density_m3 / 1000  # kg/L
# Formula for C_solution: C_soil / (Kd + theta_w / rho_b)
c_solution = c_soil / (kd + (theta_w / bulk_density_L))

# --- Step 5: Calculate PFHxS concentration in plants ---
# Assumed model: C_plant = C_solution * TSCF * PUF
c_fruit = c_solution * tscf * puf_fruit
c_legume = c_solution * tscf * puf_legume

# --- Step 6: Calculate absorbed daily intake (DI) for each food ---
# Convert intake from g to kg
intake_fruit_kg = intake_fruit_g / 1000
intake_legume_kg = intake_legume_g / 1000

# Absorbed daily intake (ug/day) = C_plant * Intake_Rate * Bioavailability
di_fruit_absorbed = c_fruit * intake_fruit_kg * bf_fruit
di_legume_absorbed = c_legume * intake_legume_kg * bf_legume

# --- Step 7: Calculate Total Dose and Hazard Quotient ---
# Total absorbed daily intake
total_di_absorbed = di_fruit_absorbed + di_legume_absorbed

# Average Daily Dose (ADD) in ug/kg/day
add = total_di_absorbed / body_weight

# Hazard Quotient (HQ)
hq = add / ref_dose

# --- Step 8: Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"1. Total PFHxS Mass Applied: {total_pfhxs_mass:,.0f} µg")
print(f"2. Total Soil Mass Affected: {total_soil_mass_kg:,.0f} kg")
print(f"3. PFHxS Concentration in Soil (C_soil): {c_soil:.4f} µg/kg")
print(f"   (Using an assumed log Koc of {log_koc_pfhxs})")
print(f"4. Soil-Water Partition Coefficient (Kd): {kd:.4f} L/kg")
print(f"5. PFHxS Concentration in Soil Solution (C_solution): {c_solution:.4f} µg/L")
print(f"6. PFHxS Concentration in Fruits: {c_fruit:.4f} µg/kg")
print(f"7. PFHxS Concentration in Legumes: {c_legume:.4f} µg/kg")
print(f"8. Absorbed Daily Intake (from fruits): {di_fruit_absorbed:.4f} µg/day")
print(f"9. Absorbed Daily Intake (from legumes): {di_legume_absorbed:.4f} µg/day")
print(f"10. Total Absorbed Dose (from all food): {total_di_absorbed:.4f} µg/day")
print("-" * 25)
print("--- Final Hazard Quotient Calculation ---")
print(f"Average Daily Dose (ADD) = Total Absorbed Dose ({total_di_absorbed:.6f} µg/day) / Body Weight ({body_weight} kg)")
print(f"ADD = {add:.6f} µg/kg/day")
print("")
print(f"Hazard Quotient (HQ) = Average Daily Dose ({add:.6f} µg/kg/day) / Reference Dose ({ref_dose} µg/kg/day)")
print(f"HQ = {hq:.6f}")

print("\nFinal Hazard Quotient:")
print(hq)
<<<0.025720>>>