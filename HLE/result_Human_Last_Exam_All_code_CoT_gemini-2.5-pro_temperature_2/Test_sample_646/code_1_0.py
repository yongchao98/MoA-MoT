import math

# --- Step 0: Define Initial Parameters ---
# Soil and Contamination Data
area = 250000  # m^2
soil_depth = 0.6  # m
f_oc = 0.03  # 3% organic carbon content
theta_w = 0.35  # L water/L soil
rho_b = 1500  # kg/m^3
foam_volume = 1000  # L
foam_conc = 1000000  # ug/L

# Human Exposure Data
bw = 80  # kg (body weight)
rfd = 0.02  # ug/kg/day (reference dose)

# Fruit Consumption Data
ir_fruit = 300  # g/day
baf_fruit = 0.5
tscf_fruit = 5.0

# Legume Consumption Data
ir_legume = 50  # g/day
baf_legume = 0.3
tscf_legume = 5.0

# Assumption for Koc (Organic Carbon-Water Partition Coefficient) for PFHxS
# This value is not given, so a typical literature value is used.
log_koc = 3.0
koc = 10**log_koc

print("--- Hazard Quotient Calculation for PFHxS Exposure ---")
print(f"Key Assumption: The organic carbon-water partition coefficient (Koc) for PFHxS is {koc} L/kg (log Koc = {log_koc}).")

# --- Step 1: Calculate Total Mass of PFHxS ---
total_pfhxs_mass = foam_volume * foam_conc
print("\nStep 1: Calculate Total Mass of PFHxS")
print(f"Total Mass = Foam Volume * Foam Concentration = {foam_volume} L * {foam_conc:,} µg/L = {total_pfhxs_mass:,.0f} µg")

# --- Step 2: Calculate Total Mass of Contaminated Soil ---
soil_volume = area * soil_depth
soil_mass = soil_volume * rho_b
print("\nStep 2: Calculate Total Mass of Contaminated Soil")
print(f"Soil Volume = Area * Depth = {area:,} m² * {soil_depth} m = {soil_volume:,.0f} m³")
print(f"Soil Mass = Soil Volume * Bulk Density = {soil_volume:,.0f} m³ * {rho_b} kg/m³ = {soil_mass:,.0f} kg")

# --- Step 3: Calculate PFHxS Concentration in Soil (C_soil) ---
c_soil = total_pfhxs_mass / soil_mass
print("\nStep 3: Calculate PFHxS Concentration in Soil (C_soil)")
print(f"C_soil = Total PFHxS Mass / Total Soil Mass = {total_pfhxs_mass:,.0f} µg / {soil_mass:,.0f} kg = {c_soil:.4f} µg/kg")

# --- Step 4: Calculate Soil-Water Partition Coefficient (Kd) ---
kd = koc * f_oc
print("\nStep 4: Calculate Soil-Water Partition Coefficient (Kd)")
print(f"Kd = Koc * f_oc = {koc} L/kg * {f_oc} = {kd:.2f} L/kg")

# --- Step 5: Calculate PFHxS Concentration in Soil Solution (Cw) ---
cw = c_soil / (kd + theta_w)
print("\nStep 5: Calculate PFHxS Concentration in Soil Solution (Cw)")
print("Using the model: Cw = C_soil / (Kd + θ_w)")
print(f"Cw = {c_soil:.4f} µg/kg / ({kd:.2f} L/kg + {theta_w} L/L) = {cw:.4f} µg/L")

# --- Step 6: Calculate PFHxS Concentration in Foods ---
# The calculation assumes plant density is 1000 g/L to convert Cw (ug/L) to C_food (ug/g)
c_fruit = (cw * tscf_fruit) / 1000
c_legume = (cw * tscf_legume) / 1000
print("\nStep 6: Calculate PFHxS Concentration in Foods (µg/g fresh weight)")
print("Using the Transpiration Stream Concentration Factor (TSCF) model: C_food = (Cw * TSCF) / 1000")
print(f"C_fruit = ({cw:.4f} µg/L * {tscf_fruit}) / 1000 = {c_fruit:.6f} µg/g")
print(f"C_legume = ({cw:.4f} µg/L * {tscf_legume}) / 1000 = {c_legume:.6f} µg/g")

# --- Step 7: Calculate Average Daily Dose (ADD) ---
absorbed_intake_fruit = c_fruit * ir_fruit * baf_fruit
absorbed_intake_legume = c_legume * ir_legume * baf_legume
total_absorbed_intake = absorbed_intake_fruit + absorbed_intake_legume
add = total_absorbed_intake / bw
print("\nStep 7: Calculate Average Daily Dose (ADD)")
print(f"Absorbed intake from fruits = {c_fruit:.6f} µg/g * {ir_fruit} g/day * {baf_fruit} = {absorbed_intake_fruit:.4f} µg/day")
print(f"Absorbed intake from legumes = {c_legume:.6f} µg/g * {ir_legume} g/day * {baf_legume} = {absorbed_intake_legume:.4f} µg/day")
print(f"Total absorbed daily intake = {absorbed_intake_fruit:.4f} + {absorbed_intake_legume:.4f} = {total_absorbed_intake:.4f} µg/day")
print(f"ADD = Total Absorbed Intake / Body Weight = {total_absorbed_intake:.4f} µg/day / {bw} kg = {add:.6f} µg/kg/day")

# --- Step 8: Calculate the Hazard Quotient (HQ) ---
hq = add / rfd
print("\nStep 8: Calculate the Hazard Quotient (HQ)")
print(f"The equation for HQ is: (Total Absorbed Daily Intake / Body Weight) / Reference Dose")
print(f"HQ = (({c_fruit:.6f} µg/g * {ir_fruit} g/day * {baf_fruit} + {c_legume:.6f} µg/g * {ir_legume} g/day * {baf_legume}) / {bw} kg) / {rfd} µg/kg/day")
print(f"HQ = {hq:.4f}")
