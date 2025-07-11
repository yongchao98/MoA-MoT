import math

# --- Given Data ---
# Soil and Contamination
area = 250000  # m^2
depth = 0.6  # m
f_oc = 0.03  # organic carbon content (3%)
theta_w = 0.35  # volumetric water content (L water/L soil)
rho_b_m3 = 1500  # bulk density (kg/m^3)
foam_volume = 1000  # L
pfhxs_conc_in_foam = 1000000  # ug/L

# Exposure Scenario
body_weight = 80  # kg
# Fruits
ir_fruit = 300  # g/day
bf_fruit = 0.5  # bioavailability factor
puf_fruit = 0.1 # plant uptake factor
# Legumes
ir_legume = 50  # g/day
bf_legume = 0.3  # bioavailability factor
puf_legume = 0.2 # plant uptake factor
# Common factors
tscf = 5 # transpiration stream concentration factor

# Toxicological Data
rfd = 0.02  # reference dose (ug/kg body weight/day)

# --- Assumptions ---
# The organic carbon-water partition coefficient (K_oc) for PFHxS is not given.
# A scientifically accepted value will be used.
k_oc = 250  # L/kg

# --- Calculations ---

print("### Hazard Quotient Calculation ###\n")

# Step 1: Calculate PFHxS concentration in soil (Cs)
print("--- Step 1: Calculate PFHxS concentration in soil (Cs) ---")
total_pfhxs = foam_volume * pfhxs_conc_in_foam
print(f"Total PFHxS added = {foam_volume} L * {pfhxs_conc_in_foam} μg/L = {total_pfhxs:,.0f} μg")

soil_volume_m3 = area * depth
soil_mass_kg = soil_volume_m3 * rho_b_m3
print(f"Total soil mass = ({area} m² * {depth} m) * {rho_b_m3} kg/m³ = {soil_mass_kg:,.0f} kg")

cs = total_pfhxs / soil_mass_kg
print(f"Concentration in soil (Cs) = {total_pfhxs:,.0f} μg / {soil_mass_kg:,.0f} kg = {cs:.4f} μg/kg\n")


# Step 2: Calculate PFHxS concentration in soil solution (Cw)
print("--- Step 2: Calculate PFHxS concentration in soil solution (Cw) ---")
print(f"Assumption: The organic carbon-water partition coefficient (K_oc) for PFHxS is {k_oc} L/kg.")
kd = k_oc * f_oc
print(f"Soil-water partition coefficient (Kd) = {k_oc} L/kg * {f_oc} = {kd:.2f} L/kg")

rho_b_kg_L = rho_b_m3 / 1000 # convert kg/m^3 to kg/L
cw = cs / (kd + (theta_w / rho_b_kg_L))
print(f"Concentration in soil water (Cw) = {cs:.4f} μg/kg / ({kd:.2f} L/kg + ({theta_w} / {rho_b_kg_L} kg/L)) = {cw:.4f} μg/L\n")


# Step 3: Calculate PFHxS concentration in produce (C_produce)
print("--- Step 3: Calculate PFHxS concentration in produce ---")
print("Method: C_produce = Cw * TSCF * PUF")
# Concentration in fruits
c_fruit = cw * tscf * puf_fruit
print(f"Concentration in fruits = {cw:.4f} μg/L * {tscf} * {puf_fruit} = {c_fruit:.4f} μg/kg (assuming 1 L ≈ 1 kg)")
# Concentration in legumes
c_legume = cw * tscf * puf_legume
print(f"Concentration in legumes = {cw:.4f} μg/L * {tscf} * {puf_legume} = {c_legume:.4f} μg/kg (assuming 1 L ≈ 1 kg)\n")


# Step 4: Calculate total daily intake (TDI)
print("--- Step 4: Calculate Total Daily Intake (TDI) ---")
# Convert plant concentration from ug/kg to ug/g
c_fruit_ugg = c_fruit / 1000
c_legume_ugg = c_legume / 1000

# Absorbed dose from each food source
di_fruit = c_fruit_ugg * ir_fruit * bf_fruit
print(f"Absorbed dose from fruits = {c_fruit_ugg:.6f} μg/g * {ir_fruit} g/day * {bf_fruit} = {di_fruit:.4f} μg/day")

di_legume = c_legume_ugg * ir_legume * bf_legume
print(f"Absorbed dose from legumes = {c_legume_ugg:.6f} μg/g * {ir_legume} g/day * {bf_legume} = {di_legume:.4f} μg/day")

total_absorbed_dose = di_fruit + di_legume
tdi = total_absorbed_dose / body_weight
print(f"Total Daily Intake (TDI) = ({di_fruit:.4f} μg/day + {di_legume:.4f} μg/day) / {body_weight} kg = {tdi:.6f} μg/kg/day\n")


# Step 5: Calculate the Hazard Quotient (HQ)
print("--- Step 5: Calculate the Hazard Quotient (HQ) ---")
hq = tdi / rfd
print(f"Hazard Quotient (HQ) = TDI / RfD")
print(f"HQ = {tdi:.6f} μg/kg/day / {rfd} μg/kg/day = {hq:.4f}\n")

print("The final Hazard Quotient is calculated to be:")
print(f"<<<{hq:.4f}>>>")