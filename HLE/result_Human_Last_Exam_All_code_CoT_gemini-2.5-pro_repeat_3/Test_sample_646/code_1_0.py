import math

# --- 1. Define Given Data ---
# Soil and Contamination Data
area = 250000  # m^2
depth = 0.6  # m
bulk_density = 1500  # kg/m^3
foc = 0.03  # fraction of organic carbon (3%)
theta_w = 0.35  # volumetric water content (L water/L soil)
foam_volume = 1000  # L
pfhxs_in_foam = 1000000  # ug/L
log_koc_pfhxs = 2.7 # Assumed literature value for PFHxS (L/kg)

# Exposure Data
body_weight = 80  # kg
rfd = 0.02  # ug/kg body weight per day

# Fruit Data
ir_fruit_g = 300  # g/day
bf_fruit = 0.5
puf_fruit = 0.1
tscf_fruit = 5

# Legume Data
ir_legume_g = 50  # g/day
bf_legume = 0.3
puf_legume = 0.2
tscf_legume = 5

print("--- Step 1: Calculate Total PFHxS and Soil Mass ---")
# Total mass of PFHxS applied
total_pfhxs_mass = foam_volume * pfhxs_in_foam
print(f"Total mass of PFHxS applied: {total_pfhxs_mass:,.0f} µg")

# Volume and mass of contaminated soil
soil_volume = area * depth
soil_mass = soil_volume * bulk_density
print(f"Total mass of contaminated soil: {soil_mass:,.0f} kg")

print("\n--- Step 2: Calculate Concentration in Soil and Soil Water (Cw) ---")
# Total concentration of PFHxS per volume of bulk soil (C_total)
c_total_ug_m3 = total_pfhxs_mass / soil_volume
# Convert to ug/L (since 1 m^3 = 1000 L)
c_total_ug_L = c_total_ug_m3 / 1000
print(f"Total PFHxS concentration (C_total) in bulk soil: {c_total_ug_L:.4f} µg/L")

# Calculate partition coefficients
koc = 10**log_koc_pfhxs
kd = koc * foc
print(f"Soil-water partition coefficient (Kd): {kd:.4f} L/kg")

# Convert bulk density to kg/L
bulk_density_kg_L = bulk_density / 1000

# Calculate concentration in soil water (Cw) using the formula: Cw = C_total / (θw + ρb * Kd)
cw = c_total_ug_L / (theta_w + bulk_density_kg_L * kd)
print(f"Concentration of PFHxS in soil water (Cw): {cw:.4f} µg/L")

print("\n--- Step 3: Calculate Concentration in Produce ---")
# Concentration in fruit (assuming 1L of plant tissue weighs 1kg)
c_fruit = cw * tscf_fruit * puf_fruit
print(f"Concentration in fruits: {c_fruit:.4f} µg/kg")

# Concentration in legumes
c_legume = cw * tscf_legume * puf_legume
print(f"Concentration in legumes: {c_legume:.4f} µg/kg")

print("\n--- Step 4: Calculate Daily Intake (DI) ---")
# Convert daily intake from g to kg
ir_fruit_kg = ir_fruit_g / 1000
ir_legume_kg = ir_legume_g / 1000

# Daily intake from fruits
di_fruit = (c_fruit * ir_fruit_kg * bf_fruit) / body_weight
print(f"Daily intake from fruits: {di_fruit:.6f} µg/kg body weight/day")

# Daily intake from legumes
di_legume = (c_legume * ir_legume_kg * bf_legume) / body_weight
print(f"Daily intake from legumes: {di_legume:.6f} µg/kg body weight/day")

print("\n--- Step 5: Calculate Total Daily Intake (TDI) and Hazard Quotient (HQ) ---")
# Total Daily Intake (TDI)
tdi = di_fruit + di_legume
print(f"Total Daily Intake (TDI): {tdi:.6f} µg/kg/day")
print(f"Reference Dose (RfD): {rfd} µg/kg/day")

# Hazard Quotient (HQ)
hq = tdi / rfd
print("\nFinal Hazard Quotient Calculation:")
print(f"HQ = TDI / RfD")
print(f"HQ = {tdi:.6f} / {rfd} = {hq:.4f}")
print("<<<0.0164>>>")