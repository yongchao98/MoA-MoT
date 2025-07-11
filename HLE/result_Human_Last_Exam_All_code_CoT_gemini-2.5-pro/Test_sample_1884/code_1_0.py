import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# --- Plan ---
# 1. Define initial parameters from the problem description.
# 2. Calculate the remaining mass of CH4 after 3 years, applying the specified oxidation rates year by year.
# 3. Convert the final mass of CH4 from metric tons to moles.
# 4. Calculate the total moles of air in the atmosphere.
# 5. Compute the increase in CH4 concentration as a mole fraction.
# 6. Convert the mole fraction to parts per billion (ppb) for the final answer.
# 7. Print all calculation steps clearly.

# --- Constants ---
total_release_ch4_tons = 250000.0  # metric tons
molar_mass_ch4_g_per_mol = 16.0
total_atmo_mass_kg = 5.1e18
# Standard value for the average molar mass of dry air
avg_molar_mass_air_g_per_mol = 28.97
KG_PER_TON = 1000.0
G_PER_KG = 1000.0

# --- Step 1: Methane mass decay over 3 years ---
print("Step 1: Calculating the final mass of methane after 3 years of oxidation.")
print(f"Initial release: {total_release_ch4_tons:,.0f} metric tons")
print("-" * 50)

# Year 1 oxidation
ch4_in_tropo_y1 = total_release_ch4_tons * 0.80
oxidation_y1_tons = ch4_in_tropo_y1 * 0.05
mass_after_y1 = total_release_ch4_tons - oxidation_y1_tons
print(f"Year 1: 5% oxidation occurs on the 80% ({ch4_in_tropo_y1:,.0f} tons) mixed into the troposphere.")
print(f"Mass oxidized in Year 1 = {ch4_in_tropo_y1:,.0f} tons * 0.05 = {oxidation_y1_tons:,.0f} tons")
print(f"Remaining mass after Year 1 = {total_release_ch4_tons:,.0f} tons - {oxidation_y1_tons:,.0f} tons = {mass_after_y1:,.0f} tons")
print("-" * 50)

# Year 2 oxidation
oxidation_y2_tons = mass_after_y1 * 0.03
mass_after_y2 = mass_after_y1 - oxidation_y2_tons
print(f"Year 2: 3% oxidation occurs on the remaining total mass.")
print(f"Mass oxidized in Year 2 = {mass_after_y1:,.0f} tons * 0.03 = {oxidation_y2_tons:,.2f} tons")
print(f"Remaining mass after Year 2 = {mass_after_y1:,.0f} tons - {oxidation_y2_tons:,.2f} tons = {mass_after_y2:,.2f} tons")
print("-" * 50)

# Year 3 oxidation
oxidation_y3_tons = mass_after_y2 * 0.03
final_ch4_mass_tons = mass_after_y2 - oxidation_y3_tons
print(f"Year 3: 3% oxidation occurs on the remaining total mass.")
print(f"Mass oxidized in Year 3 = {mass_after_y2:,.2f} tons * 0.03 = {oxidation_y3_tons:,.2f} tons")
print(f"Final remaining mass after 3 years = {mass_after_y2:,.2f} tons - {oxidation_y3_tons:,.2f} tons = {final_ch4_mass_tons:,.2f} tons")
print("-" * 50)

# --- Step 2 & 3: Convert to moles and calculate atmospheric moles ---
print("\nStep 2: Calculating moles of methane and atmosphere.")
# Moles of CH4
final_ch4_mass_g = final_ch4_mass_tons * KG_PER_TON * G_PER_KG
moles_ch4 = final_ch4_mass_g / molar_mass_ch4_g_per_mol
print(f"Final CH4 mass in grams = {final_ch4_mass_tons:,.2f} tons * {KG_PER_TON} kg/ton * {G_PER_KG} g/kg = {final_ch4_mass_g:.3e} g")
print(f"Moles of CH4 = {final_ch4_mass_g:.3e} g / {molar_mass_ch4_g_per_mol} g/mol = {moles_ch4:.3e} moles")
print("-" * 50)

# Moles of Atmosphere
total_atmo_mass_g = total_atmo_mass_kg * G_PER_KG
total_moles_atmo = total_atmo_mass_g / avg_molar_mass_air_g_per_mol
print(f"Total atmosphere mass in grams = {total_atmo_mass_kg:.1e} kg * {G_PER_KG} g/kg = {total_atmo_mass_g:.1e} g")
print(f"Total moles of atmosphere = {total_atmo_mass_g:.1e} g / {avg_molar_mass_air_g_per_mol} g/mol = {total_moles_atmo:.3e} moles")
print("-" * 50)

# --- Step 4: Calculate ppb increase ---
mole_fraction_increase = moles_ch4 / total_moles_atmo
ppb_increase = mole_fraction_increase * 1e9

print("\nStep 3: Final Calculation for Concentration Increase (ppb).")
print(f"Mole fraction increase = Moles CH4 / Total Moles Atmosphere")
print(f"Mole fraction = {moles_ch4:.3e} moles / {total_moles_atmo:.3e} moles = {mole_fraction_increase:.3e}")
print(f"Concentration increase (ppb) = Mole Fraction * 1,000,000,000")
print(f"Final ppb increase = {mole_fraction_increase:.3e} * 1e9 = {ppb_increase:.3f} ppb")
print("-" * 50)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = buffer.getvalue()
# Print the output to the user
print(output)

print(f"<<<{ppb_increase:.3f}>>>")