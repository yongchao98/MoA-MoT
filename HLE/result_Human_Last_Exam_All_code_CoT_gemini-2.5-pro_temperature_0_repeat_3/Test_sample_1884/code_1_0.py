import sys

# Step 1: Define constants and initial values
total_ch4_leak_metric_tons = 250000
# Convert metric tons to kilograms
total_ch4_leak_kg = total_ch4_leak_metric_tons * 1000
# Molecular weight of CH4 in kg/mol
mw_ch4_kg_mol = 16 / 1000
# Total mass of the atmosphere in kg
mass_atmosphere_kg = 5.1 * 10**18
# Average molecular weight of air in kg/mol (approx. 29 g/mol)
mw_air_kg_mol = 29 / 1000

# Step 2: Calculate the initial total moles of CH4 released
initial_moles_ch4 = total_ch4_leak_kg / mw_ch4_kg_mol
print(f"Initial total moles of CH4 released: {initial_moles_ch4:.4e} mol")

# Step 3: Distribute methane into troposphere and stratosphere
moles_ch4_tropo = initial_moles_ch4 * 0.80
moles_ch4_strato = initial_moles_ch4 * 0.20
print(f"Initial CH4 in troposphere: {moles_ch4_tropo:.4e} mol")
print(f"Initial CH4 in stratosphere: {moles_ch4_strato:.4e} mol\n")

# Step 4: Apply annual oxidation reductions sequentially
# Year 1: 5% reduction in troposphere, 0% in stratosphere
moles_ch4_tropo_yr1 = moles_ch4_tropo * (1 - 0.05)
moles_ch4_strato_yr1 = moles_ch4_strato # No change in year 1
print("--- After Year 1 ---")
print(f"Remaining CH4 in troposphere: {moles_ch4_tropo_yr1:.4e} mol")
print(f"Remaining CH4 in stratosphere: {moles_ch4_strato_yr1:.4e} mol\n")

# Year 2: 3% reduction in both layers
moles_ch4_tropo_yr2 = moles_ch4_tropo_yr1 * (1 - 0.03)
moles_ch4_strato_yr2 = moles_ch4_strato_yr1 * (1 - 0.03)
print("--- After Year 2 ---")
print(f"Remaining CH4 in troposphere: {moles_ch4_tropo_yr2:.4e} mol")
print(f"Remaining CH4 in stratosphere: {moles_ch4_strato_yr2:.4e} mol\n")

# Year 3: 3% reduction in both layers
moles_ch4_tropo_yr3 = moles_ch4_tropo_yr2 * (1 - 0.03)
moles_ch4_strato_yr3 = moles_ch4_strato_yr2 * (1 - 0.03)
print("--- After Year 3 ---")
print(f"Remaining CH4 in troposphere: {moles_ch4_tropo_yr3:.4e} mol")
print(f"Remaining CH4 in stratosphere: {moles_ch4_strato_yr3:.4e} mol\n")

# Step 5: Calculate the total final moles of CH4
final_moles_ch4 = moles_ch4_tropo_yr3 + moles_ch4_strato_yr3

# Step 6: Calculate the total moles of air in the atmosphere
total_moles_air = mass_atmosphere_kg / mw_air_kg_mol

# Step 7: Calculate the final concentration increase in ppb
# ppb = (moles_of_gas / total_moles_of_air) * 10^9
final_ppb_increase = (final_moles_ch4 / total_moles_air) * 10**9

print("--- Final Calculation ---")
print(f"Final total moles of CH4 after 3 years: {final_moles_ch4:.4e} mol")
print(f"Total moles of air in the atmosphere: {total_moles_air:.4e} mol")
print("\nFinal Equation: (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
print(f"({final_moles_ch4:.4e} / {total_moles_air:.4e}) * 1e9")
print(f"\nCalculated increase in CH4 concentration after 3 years: {final_ppb_increase:.4f} ppb")

# Final answer block
sys.stdout.write(f"<<<{final_ppb_increase:.4f}>>>")