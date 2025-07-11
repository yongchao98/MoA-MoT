import math

# --- Given Data ---

# Soil Properties
area_m2 = 250000
depth_m = 0.6
f_oc = 0.03  # organic carbon content fraction
theta_w = 0.35  # volumetric water content (L water/L soil)
rho_b_kg_m3 = 1500  # bulk density (kg/m^3)

# Contaminant Properties
foam_volume_L = 1000
pfhxs_conc_in_foam_ug_L = 1000000
# Assumption: The organic carbon-water partitioning coefficient (Koc) for PFHxS is not given.
# A literature value of log Koc = 2.8 is used.
log_koc = 2.8
k_oc_L_kg = 10**log_koc

# Human Exposure Factors
body_weight_kg = 80
rfd_ug_kg_day = 0.02  # Reference Dose

# Fruit Consumption
fruit_intake_g_day = 300
fruit_bf = 0.5  # bioavailability factor
fruit_puf = 0.1  # plant uptake factor
fruit_tscf = 5  # transpiration stream concentration factor

# Legume Consumption
legume_intake_g_day = 50
legume_bf = 0.3  # bioavailability factor
legume_puf = 0.2  # plant uptake factor
legume_tscf = 5  # transpiration stream concentration factor

# --- Calculations ---

print("Step-by-step Calculation of the Hazard Quotient (HQ)\n")

# 1. Total mass of PFHxS applied
total_pfhxs_mass_ug = foam_volume_L * pfhxs_conc_in_foam_ug_L
print(f"1. Total mass of PFHxS applied = {foam_volume_L} L * {pfhxs_conc_in_foam_ug_L} μg/L = {total_pfhxs_mass_ug:e} μg")

# 2. Total mass of contaminated soil
soil_volume_m3 = area_m2 * depth_m
total_soil_mass_kg = soil_volume_m3 * rho_b_kg_m3
print(f"2. Total mass of soil = ({area_m2} m² * {depth_m} m) * {rho_b_kg_m3} kg/m³ = {total_soil_mass_kg:e} kg")

# 3. Concentration of PFHxS in soil (C_soil)
c_soil_ug_kg = total_pfhxs_mass_ug / total_soil_mass_kg
print(f"3. Concentration in soil (C_soil) = {total_pfhxs_mass_ug:e} μg / {total_soil_mass_kg:e} kg = {c_soil_ug_kg:.4f} μg/kg")

# 4. Concentration of PFHxS in soil solution (C_solution)
# 4a. Calculate soil-water partition coefficient (Kd)
k_d_L_kg = k_oc_L_kg * f_oc
print(f"4a. Soil-water partition coefficient (Kd) = {k_oc_L_kg:.2f} L/kg * {f_oc} = {k_d_L_kg:.4f} L/kg")

# 4b. Convert bulk density to kg/L for unit consistency
rho_b_kg_L = rho_b_kg_m3 / 1000
# 4c. Calculate C_solution
c_solution_ug_L = c_soil_ug_kg / (k_d_L_kg + (theta_w / rho_b_kg_L))
print(f"4b. Concentration in soil solution (C_solution) = {c_soil_ug_kg:.4f} μg/kg / ({k_d_L_kg:.4f} L/kg + ({theta_w} / {rho_b_kg_L} kg/L)) = {c_solution_ug_L:.4f} μg/L")

# 5. Absorbed Daily Dose (ADD) from Fruits
# 5a. Concentration in fruits (assuming plant density ~ 1 kg/L)
c_plant_fruit_ug_kg = c_solution_ug_L * fruit_tscf * fruit_puf
print(f"\n--- Fruit Intake ---")
print(f"5a. Concentration in fruits (C_fruit) = {c_solution_ug_L:.4f} μg/L * {fruit_tscf} * {fruit_puf} = {c_plant_fruit_ug_kg:.4f} μg/kg")
# 5b. Daily intake rate in kg/day
fruit_intake_kg_day = fruit_intake_g_day / 1000
# 5c. ADD from fruits
add_fruit_ug_kg_day = (c_plant_fruit_ug_kg * fruit_intake_kg_day * fruit_bf) / body_weight_kg
print(f"5b. Absorbed Daily Dose (ADD_fruit) = ({c_plant_fruit_ug_kg:.4f} μg/kg * {fruit_intake_kg_day} kg/day * {fruit_bf}) / {body_weight_kg} kg = {add_fruit_ug_kg_day:.6f} μg/kg/day")

# 6. Absorbed Daily Dose (ADD) from Legumes
# 6a. Concentration in legumes (assuming plant density ~ 1 kg/L)
c_plant_legume_ug_kg = c_solution_ug_L * legume_tscf * legume_puf
print(f"\n--- Legume Intake ---")
print(f"6a. Concentration in legumes (C_legume) = {c_solution_ug_L:.4f} μg/L * {legume_tscf} * {legume_puf} = {c_plant_legume_ug_kg:.4f} μg/kg")
# 6b. Daily intake rate in kg/day
legume_intake_kg_day = legume_intake_g_day / 1000
# 6c. ADD from legumes
add_legume_ug_kg_day = (c_plant_legume_ug_kg * legume_intake_kg_day * legume_bf) / body_weight_kg
print(f"6b. Absorbed Daily Dose (ADD_legume) = ({c_plant_legume_ug_kg:.4f} μg/kg * {legume_intake_kg_day} kg/day * {legume_bf}) / {body_weight_kg} kg = {add_legume_ug_kg_day:.6f} μg/kg/day")

# 7. Total Absorbed Daily Dose
total_add_ug_kg_day = add_fruit_ug_kg_day + add_legume_ug_kg_day
print(f"\n--- Total Exposure and Hazard ---")
print(f"7. Total ADD = {add_fruit_ug_kg_day:.6f} μg/kg/day + {add_legume_ug_kg_day:.6f} μg/kg/day = {total_add_ug_kg_day:.6f} μg/kg/day")

# 8. Hazard Quotient (HQ)
hq = total_add_ug_kg_day / rfd_ug_kg_day
print(f"\n8. Final Hazard Quotient (HQ) = Total ADD / Reference Dose (RfD)")
print(f"   HQ = {total_add_ug_kg_day:.6f} / {rfd_ug_kg_day} = {hq:.4f}")

print(f"\nThe final Hazard Quotient is {hq:.4f}.")
print("<<<" + str(round(hq, 4)) + ">>>")