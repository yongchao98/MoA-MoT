import math

# Step 1: Define initial parameters from the problem description
# --- Contaminant details ---
foam_volume_L = 1000  # L
pfhxs_in_foam_ug_per_L = 1000000  # µg/L

# --- Soil characteristics ---
area_m2 = 250000  # m²
depth_m = 0.6  # m
foc = 0.03  # 3% organic carbon content
theta_w = 0.35  # L water/L soil, volumetric water content
rho_b_kg_m3 = 1500  # kg/m³

# --- Human exposure details ---
body_weight_kg = 80  # kg
rfd_ug_per_kg_day = 0.02  # Reference Dose in µg/kg/day

# --- Fruit consumption details ---
ir_fruit_g_day = 300  # g/day
b_fruit = 0.5  # bioavailability factor
tscf_fruit = 5 # transpiration stream concentration factor

# --- Legume consumption details ---
ir_legume_g_day = 50  # g/day
b_legume = 0.3  # bioavailability factor
tscf_legume = 5 # transpiration stream concentration factor

# --- Assumed parameter ---
# Koc for PFHxS is not given, a literature value of log(Koc)=2.5 is used.
koc_L_per_kg = 10**2.5

# Step 2: Calculate total mass of PFHxS applied
total_pfhxs_mass_ug = foam_volume_L * pfhxs_in_foam_ug_per_L

# Step 3: Calculate volume and mass of contaminated soil
soil_volume_m3 = area_m2 * depth_m
soil_mass_kg = soil_volume_m3 * rho_b_kg_m3

# Step 4: Calculate PFHxS concentration in total soil (mass/mass)
cs_ug_per_kg = total_pfhxs_mass_ug / soil_mass_kg

# Step 5: Calculate PFHxS concentration in soil pore water (Cw)
# Calculate soil partition coefficient (Kd)
kd_L_per_kg = koc_L_per_kg * foc

# Convert bulk density to kg/L for unit consistency
rho_b_kg_L = rho_b_kg_m3 / 1000

# Calculate Cw using the partitioning equation
# Cw = Cs / (Kd + θw/ρb)
cw_ug_per_L = cs_ug_per_kg / (kd_L_per_kg + (theta_w / rho_b_kg_L))

# Step 6: Calculate PFHxS concentration in plants
# Since TSCF is the same for both, plant concentrations are the same.
c_plant_ug_per_L = cw_ug_per_L * tscf_fruit # Using fruit TSCF, but it's the same as legume
# Convert plant concentration to mass basis (µg/g), assuming plant density is 1 g/mL or 1000 g/L
c_plant_ug_per_g = c_plant_ug_per_L / 1000

# Step 7: Calculate daily intake (DI) for each food type
di_fruit_ug_kg_day = (c_plant_ug_per_g * ir_fruit_g_day * b_fruit) / body_weight_kg
di_legume_ug_kg_day = (c_plant_ug_per_g * ir_legume_g_day * b_legume) / body_weight_kg

# Calculate Total Daily Intake (TDI)
tdi_ug_kg_day = di_fruit_ug_kg_day + di_legume_ug_kg_day

# Step 8: Calculate the Hazard Quotient (HQ)
hq = tdi_ug_kg_day / rfd_ug_per_kg_day

# Final Output: Print the calculation of the final equation as requested
print("The hazard quotient (HQ) is calculated as follows:")
print("HQ = (Daily Intake from Fruit + Daily Intake from Legumes) / Reference Dose")
print(f"HQ = ({di_fruit_ug_kg_day:.6f} + {di_legume_ug_kg_day:.6f}) / {rfd_ug_per_kg_day}")
print(f"HQ = {tdi_ug_kg_day:.6f} / {rfd_ug_per_kg_day}")
print(f"Hazard Quotient (HQ) = {hq:.4f}")

# The final answer in the requested format
print(f"\n<<<{hq:.4f}>>>")