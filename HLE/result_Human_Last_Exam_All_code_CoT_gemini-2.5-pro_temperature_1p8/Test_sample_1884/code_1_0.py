import sys

# Step 1: Define the given constants. A key constant, the average molar mass of air 
# (approx. 29 g/mol), is not given but is required to convert from mass to mole fraction (ppb). 
# We will use this standard value.
total_methane_release_tons = 250000
total_methane_release_kg = total_methane_release_tons * 1000  # Convert metric tons to kg

mass_atmosphere_kg = 5.1e18  # Total mass of the atmosphere
mm_ch4_g_mol = 16.0
mm_ch4_kg_mol = mm_ch4_g_mol / 1000  # Molar mass of CH4 in kg/mol

# Standard average molar mass of dry air
mm_air_g_mol = 29.0
mm_air_kg_mol = mm_air_g_mol / 1000  # Molar mass of air in kg/mol

# Step 2: Calculate methane mass reduction over 3 years

# At the time of release, methane is split into two portions.
methane_in_troposphere_initial_kg = total_methane_release_kg * 0.80
methane_for_stratosphere_kg = total_methane_release_kg * 0.20

# During Year 1, only the tropospheric portion is oxidized.
oxidation_y1 = 0.05
methane_in_troposphere_after_y1 = methane_in_troposphere_initial_kg * (1 - oxidation_y1)

# At the end of Year 1, the total remaining methane is the sum of the oxidized tropospheric part 
# and the untouched other part.
total_methane_end_y1 = methane_in_troposphere_after_y1 + methane_for_stratosphere_kg

# During Year 2, the *total* remaining methane is oxidized by 3%.
oxidation_y2_y3 = 0.03
total_methane_end_y2 = total_methane_end_y1 * (1 - oxidation_y2_y3)

# During Year 3, the remaining amount is oxidized by 3% again.
final_methane_mass_kg = total_methane_end_y2 * (1 - oxidation_y2_y3)

# Step 3: Calculate the final concentration increase in ppb

# Convert final methane mass to moles
final_moles_ch4 = final_methane_mass_kg / mm_ch4_kg_mol

# Convert total atmospheric mass to moles
total_moles_air = mass_atmosphere_kg / mm_air_kg_mol

# Calculate the concentration increase in ppb (parts per billion)
# ppb is a mole fraction, so (moles_ch4 / total_moles_air) * 1e9
ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9

# Step 4: Print the calculation steps and the final result

print("Calculation of Methane Concentration Increase After 3 Years")
print("----------------------------------------------------------")
print(f"1. Initial methane release: {total_methane_release_kg:e} kg")
print(f"   - Amount mixed into troposphere (80%): {methane_in_troposphere_initial_kg:.0f} kg")
print(f"   - Amount set aside for stratosphere (20%): {methane_for_stratosphere_kg:.0f} kg")

print("\n2. Mass remaining after Year 1 (5% oxidation in troposphere):")
print(f"   - Total remaining CH4 at end of Year 1: ({methane_in_troposphere_initial_kg:.0f} * (1 - {oxidation_y1})) + {methane_for_stratosphere_kg:.0f} = {total_methane_end_y1:.0f} kg")

print("\n3. Mass remaining after Year 2 (3% oxidation of total):")
print(f"   - Total remaining CH4 at end of Year 2: {total_methane_end_y1:.0f} * (1 - {oxidation_y2_y3}) = {total_methane_end_y2:.0f} kg")

print("\n4. Mass remaining after Year 3 (3% oxidation of total):")
print(f"   - Final remaining CH4 mass: {total_methane_end_y2:.0f} * (1 - {oxidation_y2_y3}) = {final_methane_mass_kg:.0f} kg")

print("\n5. Final Concentration Calculation (ppb):")
print(f"   - Final moles of CH4 = {final_methane_mass_kg:.0f} kg / {mm_ch4_kg_mol} kg/mol = {final_moles_ch4:e} mol")
print(f"   - Total moles of Air = {mass_atmosphere_kg:e} kg / {mm_air_kg_mol} kg/mol = {total_moles_air:e} mol")

print("\n   Final Equation:")
print(f"   ppb Increase = (Moles CH4 / Total Moles Air) * 1,000,000,000")
print(f"   ppb Increase = ({final_moles_ch4:.4e} mol / {total_moles_air:.4e} mol) * 1e9")
print("----------------------------------------------------------")
print(f"Final calculated increase in atmospheric methane concentration: {ppb_increase:.4f} ppb")
<<<0.0803>>>