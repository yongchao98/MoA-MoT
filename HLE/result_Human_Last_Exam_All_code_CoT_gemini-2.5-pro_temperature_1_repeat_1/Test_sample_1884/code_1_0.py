import math

# Plan:
# 1. Define the given constants.
# 2. Calculate the total initial moles of methane (CH4) released.
# 3. Model the reduction of methane over 3 years based on the specified mixing and oxidation rules.
# 4. Calculate the total moles of air in the atmosphere.
# 5. Compute the final concentration increase in parts per billion (ppb).
# 6. Print the steps and the final equation.

# Step 1: Define constants
ch4_release_tons = 250000.0
ch4_molar_mass_g_mol = 16.0
atm_total_mass_kg = 5.1e18
# The average molar mass of Earth's air is approximately 29.0 g/mol.
air_molar_mass_g_mol = 29.0
tropo_mix_fraction = 0.80
strato_mix_fraction = 1.0 - tropo_mix_fraction
y1_oxidation_rate = 0.05
subsequent_oxidation_rate = 0.03
ppb_conversion_factor = 1e9

# Step 2: Calculate total initial moles of CH4 released
print("--- Initial Methane Calculation ---")
ch4_release_g = ch4_release_tons * 1e6
initial_ch4_moles = ch4_release_g / ch4_molar_mass_g_mol
print(f"Total CH4 released: {ch4_release_tons} metric tons")
print(f"Initial moles of CH4 = {ch4_release_g:.2e} g / {ch4_molar_mass_g_mol} g/mol = {initial_ch4_moles:.4e} moles\n")

# Step 3: Track methane reduction over 3 years
print("--- Methane Reduction Over 3 Years ---")
# At the time of the leak (Year 0)
ch4_moles_for_tropo = initial_ch4_moles * tropo_mix_fraction
ch4_moles_for_strato = initial_ch4_moles * strato_mix_fraction

# End of Year 1
# 5% of the tropospheric portion is oxidized.
ch4_tropo_after_y1 = ch4_moles_for_tropo * (1 - y1_oxidation_rate)
total_ch4_after_y1 = ch4_tropo_after_y1 + ch4_moles_for_strato
print(f"Year 1: Methane reduced by {y1_oxidation_rate*100}% in the troposphere.")
print(f"Remaining CH4 after Year 1 = ({ch4_moles_for_tropo:.4e} * (1 - {y1_oxidation_rate})) + {ch4_moles_for_strato:.4e} = {total_ch4_after_y1:.4e} moles")

# End of Year 2
# 3% of the total remaining methane is oxidized.
total_ch4_after_y2 = total_ch4_after_y1 * (1 - subsequent_oxidation_rate)
print(f"Year 2: Total remaining methane reduced by {subsequent_oxidation_rate*100}%.")
print(f"Remaining CH4 after Year 2 = {total_ch4_after_y1:.4e} * (1 - {subsequent_oxidation_rate}) = {total_ch4_after_y2:.4e} moles")

# End of Year 3
# Another 3% of the total remaining methane is oxidized.
final_ch4_moles = total_ch4_after_y2 * (1 - subsequent_oxidation_rate)
print(f"Year 3: Total remaining methane reduced by another {subsequent_oxidation_rate*100}%.")
print(f"Final remaining CH4 after Year 3 = {total_ch4_after_y2:.4e} * (1 - {subsequent_oxidation_rate}) = {final_ch4_moles:.4e} moles\n")

# Step 4: Calculate total moles of air
print("--- Atmospheric Calculation ---")
atm_total_mass_g = atm_total_mass_kg * 1e3
total_air_moles = atm_total_mass_g / air_molar_mass_g_mol
print(f"Total mass of atmosphere: {atm_total_mass_g:.2e} g")
print(f"Total moles of air = {atm_total_mass_g:.2e} g / {air_molar_mass_g_mol} g/mol = {total_air_moles:.4e} moles\n")

# Step 5: Compute final concentration increase in ppb
print("--- Final Concentration Calculation ---")
final_increase_ppb = (final_ch4_moles / total_air_moles) * ppb_conversion_factor
print("Final Increase (ppb) = (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
print(f"Final Increase (ppb) = ({final_ch4_moles:.5f} / {total_air_moles:.5f}) * {int(ppb_conversion_factor)}")
print(f"Result: {final_increase_ppb:.4f} ppb")

print(f"<<<{final_increase_ppb}>>>")