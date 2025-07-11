import sys

# --- Given Data ---

# Soil and Site characteristics
area = 250000  # m^2
depth = 0.6  # m
foc = 0.03  # fraction of organic carbon
theta_w = 0.35  # volumetric water content (L water/L soil)
rho_b_m3 = 1500  # bulk density (kg/m^3)

# Contamination event
foam_volume = 1000  # L
pfhxs_in_foam_conc = 1000000  # ug/L

# Human exposure factors
body_weight = 80  # kg

# Fruits consumption
fruit_intake_rate = 0.3  # kg/day (300g)
fruit_bioavailability = 0.5
fruit_puf = 0.1  # plant uptake factor
fruit_tscf = 5   # transpiration stream concentration factor

# Legumes consumption
legume_intake_rate = 0.05  # kg/day (50g)
legume_bioavailability = 0.3
legume_puf = 0.2  # plant uptake factor
legume_tscf = 5   # transpiration stream concentration factor

# Toxicological data
rfd = 0.02  # reference dose (ug/kg body weight per day)

# Assumed value (not provided in problem)
# Koc (Organic Carbon-Water Partition Coefficient) for PFHxS.
# A literature value of 1000 L/kg (log Koc = 3) is a reasonable assumption.
koc = 1000 # L/kg


# --- Calculations ---

# Step 1: Calculate total mass of PFHxS applied
total_pfhxs_mass = foam_volume * pfhxs_in_foam_conc  # in ug

# Step 2: Calculate concentration of PFHxS in soil (Cs)
soil_volume = area * depth  # in m^3
total_soil_mass = soil_volume * rho_b_m3  # in kg
cs = total_pfhxs_mass / total_soil_mass  # in ug/kg

# Step 3: Calculate concentration of PFHxS in soil solution (Cw)
kd = foc * koc  # soil-water partition coefficient (L/kg)
# Convert bulk density from kg/m^3 to kg/L for unit consistency
rho_b_L = rho_b_m3 / 1000.0  # kg/L
# Cw = Cs / (Kd + (theta_w / rho_b))
cw = cs / (kd + (theta_w / rho_b_L))  # in ug/L

# Step 4: Calculate PFHxS concentration in plants
# C_plant = Cw * TSCF * PUF. Assuming this results in ug/kg of plant tissue.
c_fruit = cw * fruit_tscf * fruit_puf  # in ug/kg
c_legume = cw * legume_tscf * legume_puf  # in ug/kg

# Step 5: Calculate total daily intake (DI) and Chronic Daily Intake (CDI)
# Daily intake from fruits, considering bioavailability
di_fruit = c_fruit * fruit_intake_rate * fruit_bioavailability  # in ug/day

# Daily intake from legumes, considering bioavailability
di_legume = c_legume * legume_intake_rate * legume_bioavailability # in ug/day

# Total daily intake
total_di = di_fruit + di_legume  # in ug/day

# Chronic Daily Intake (CDI) normalized by body weight
cdi = total_di / body_weight  # in ug/kg/day

# Step 6: Calculate Hazard Quotient (HQ)
hq = cdi / rfd

# --- Final Output ---
print("This script calculates the Hazard Quotient (HQ) for PFHxS exposure.")
print("\n--- Final Equation ---")
print(f"Hazard Quotient = Chronic Daily Intake / Reference Dose")
print(f"Hazard Quotient = {cdi:.8f} µg/kg/day / {rfd} µg/kg/day")
print(f"\nThe calculated Hazard Quotient (HQ) is: {hq:.4f}")

# Suppress the prompt ending with <<< at the beginning of a new line.
# Flush stdout to ensure the output is written before the final answer tag.
sys.stdout.flush()

print(f"<<<{hq:.4f}>>>")