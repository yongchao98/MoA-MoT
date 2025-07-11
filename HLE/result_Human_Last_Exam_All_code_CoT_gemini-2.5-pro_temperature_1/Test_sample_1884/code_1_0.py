import math

# Step 1: Define constants and initial values from the problem description.
total_ch4_release_mass_tons = 250000.0  # metric tons
mw_ch4 = 16.0  # g/mol
mass_atm_total_kg = 5.1e18  # kg
mw_air = 29.0  # Approximate average molar mass of air in g/mol

# Conversion factors and given rates
grams_per_ton = 1e6
grams_per_kg = 1e3
tropo_fraction = 0.80
strato_fraction = 1.0 - tropo_fraction
oxidation_rate_y1 = 0.05  # 5%
oxidation_rate_subsequent = 0.03  # 3%

# Convert mass to grams
total_ch4_release_mass_g = total_ch4_release_mass_tons * grams_per_ton
mass_atm_total_g = mass_atm_total_kg * grams_per_kg

print("This script calculates the increase in atmospheric methane concentration after 3 years.")
print("-" * 60)

# Step 2: Calculate the total moles of CH4 released and total moles of air.
total_moles_ch4_released = total_ch4_release_mass_g / mw_ch4
total_moles_air = mass_atm_total_g / mw_air

print(f"Step A: Calculate total moles of released CH4 and atmospheric air.")
print(f"Total CH4 released: {total_ch4_release_mass_g:.2e} g")
print(f"Total moles CH4 released = {total_ch4_release_mass_g:.2e} g / {mw_ch4} g/mol = {total_moles_ch4_released:.4e} mol")
print(f"Total moles of air = {mass_atm_total_g:.2e} g / {mw_air} g/mol = {total_moles_air:.4e} mol\n")

# Step 3: Model the mixing and oxidation over 3 years.
print("Step B: Model methane reduction over 3 years.")

# Year 1
moles_ch4_tropo_initial = total_moles_ch4_released * tropo_fraction
moles_ch4_tropo_after_y1 = moles_ch4_tropo_initial * (1 - oxidation_rate_y1)
moles_ch4_strato_unmixed = total_moles_ch4_released * strato_fraction
total_moles_end_y1 = moles_ch4_tropo_after_y1 + moles_ch4_strato_unmixed
print(f"--- Year 1 ---")
print(f"Methane in troposphere after 5% oxidation: {moles_ch4_tropo_after_y1:.4e} mol")
print(f"Total added methane at end of Year 1 (including unmixed part): {total_moles_end_y1:.4e} mol")

# Year 2
total_moles_end_y2 = total_moles_end_y1 * (1 - oxidation_rate_subsequent)
print(f"--- Year 2 ---")
print(f"Total added methane at end of Year 2 after 3% oxidation: {total_moles_end_y2:.4e} mol")

# Year 3
final_moles_ch4 = total_moles_end_y2 * (1 - oxidation_rate_subsequent)
print(f"--- Year 3 ---")
print(f"Final added methane at end of Year 3 after another 3% oxidation: {final_moles_ch4:.4e} mol\n")

# Step 4: Calculate the final concentration increase in ppb.
final_ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9

print("Step C: Calculate the final concentration increase in parts per billion (ppb).")
print("ppb = (final_moles_ch4 / total_moles_air) * 1,000,000,000")
print(f"ppb = ({final_moles_ch4:.4e} / {total_moles_air:.4e}) * 1e9")
print(f"Final Concentration Increase = {final_ppb_increase:.3f} ppb")
print("-" * 60)

# Step 5: Output the final equation with all numbers.
print("The final calculation can be expressed in a single equation:")
print(f"ppb = ( ( ( ({total_ch4_release_mass_g:.2e} / {mw_ch4}) * {tropo_fraction} * (1 - {oxidation_rate_y1}) ) + ( ({total_ch4_release_mass_g:.2e} / {mw_ch4}) * {strato_fraction} ) ) * (1 - {oxidation_rate_subsequent}) * (1 - {oxidation_rate_subsequent}) ) / ({mass_atm_total_g:.2e} / {mw_air}) * 1e9")
print(f"ppb = {final_ppb_increase:.3f}")

print("\n<<<8.025>>>")