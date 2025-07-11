import math

# Step 1: Define constants and calculate initial moles of CH4 released.
total_ch4_release_tons = 250000.0
mw_ch4_g_mol = 16.0
mass_atmosphere_kg = 5.1e18
mw_air_g_mol = 29.0

# Convert release from tons to grams
total_ch4_release_g = total_ch4_release_tons * 1000 * 1000

# Calculate total moles of CH4 released
moles_ch4_released = total_ch4_release_g / mw_ch4_g_mol

print("Step 1: Calculate Total Moles of CH4 Released")
print(f"({total_ch4_release_tons:.0f} tons * 1000 kg/ton * 1000 g/kg) / {mw_ch4_g_mol} g/mol = {moles_ch4_released:e} moles of CH4")
print("-" * 50)

# Step 2: Account for distribution and oxidation over 3 years.
# Year 1: 80% is in the troposphere and is oxidized by 5%. The other 20% is separate.
moles_in_troposphere_Y1 = moles_ch4_released * 0.80
moles_for_stratosphere = moles_ch4_released * 0.20
oxidation_Y1 = 0.05
moles_after_Y1_oxidation = moles_in_troposphere_Y1 * (1 - oxidation_Y1)
total_moles_after_Y1 = moles_after_Y1_oxidation + moles_for_stratosphere

# Year 2 and 3: The total remaining amount is oxidized by 3% each year.
oxidation_subsequent_years = 0.03
total_moles_after_Y2 = total_moles_after_Y1 * (1 - oxidation_subsequent_years)
final_moles_ch4_remaining = total_moles_after_Y2 * (1 - oxidation_subsequent_years)

print("Step 2: Calculate Moles of CH4 Remaining After 3 Years of Oxidation")
print(f"Initial moles for troposphere (80%): {moles_in_troposphere_Y1:e}")
print(f"Moles remaining after 5% oxidation in Year 1: {moles_in_troposphere_Y1:e} * (1 - {oxidation_Y1}) = {moles_after_Y1_oxidation:e}")
print(f"Total moles after Year 1 (adding back the 20%): {moles_after_Y1_oxidation:e} + {moles_for_stratosphere:e} = {total_moles_after_Y1:e}")
print(f"Total moles after Year 2 (3% oxidation): {total_moles_after_Y1:e} * (1 - {oxidation_subsequent_years}) = {total_moles_after_Y2:e}")
print(f"Final moles remaining after Year 3 (another 3% oxidation): {total_moles_after_Y2:e} * (1 - {oxidation_subsequent_years}) = {final_moles_ch4_remaining:e}")
print("-" * 50)

# Step 3: Calculate total moles of air.
mass_atmosphere_g = mass_atmosphere_kg * 1000
total_moles_air = mass_atmosphere_g / mw_air_g_mol

print("Step 3: Calculate Total Moles of Air in the Atmosphere")
print(f"({mass_atmosphere_kg:e} kg * 1000 g/kg) / {mw_air_g_mol} g/mol = {total_moles_air:e} moles of air")
print("-" * 50)

# Step 4: Calculate the final concentration increase in ppb.
final_ppb_increase = (final_moles_ch4_remaining / total_moles_air) * 1e9

print("Step 4: Calculate the Final Increase in CH4 Concentration (ppb)")
print(f"Final ppb = (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
print(f"({final_moles_ch4_remaining:e} moles CH4 / {total_moles_air:e} moles air) * 1e9 = {final_ppb_increase:.5f} ppb")
print("-" * 50)
print(f"The final calculated increase in atmospheric methane concentration after 3 years is {final_ppb_increase:.5f} ppb.")

print(f"<<<{final_ppb_increase:.5f}>>>")