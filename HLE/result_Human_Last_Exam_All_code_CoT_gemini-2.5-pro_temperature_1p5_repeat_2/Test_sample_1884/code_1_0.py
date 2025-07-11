import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string stream
sys.stdout = io.StringIO()

# --- 1. Define constants ---
# Methane leak properties
total_ch4_release_tons = 250000.0
mw_ch4 = 16.0  # g/mol

# Atmospheric properties
total_mass_atmosphere_kg = 5.1e18
mw_air = 28.97  # g/mol (standard average molecular weight of dry air)

# Mixing and Oxidation parameters
leak_to_tropo_fraction = 0.80
leak_to_strato_fraction = 1 - leak_to_tropo_fraction
years_for_strato_mix = 3.0
tropo_oxidation_rate_y1 = 0.05
annual_oxidation_rate_subsequent = 0.03

# --- 2. Calculate initial moles of CH4 released ---
# Convert metric tons to grams
total_ch4_release_g = total_ch4_release_tons * 1e6
# Calculate total moles of leaked methane
total_leaked_moles = total_ch4_release_g / mw_ch4

# --- 3. Calculate distribution of leaked CH4 per year ---
moles_mixed_to_tropo_y1 = total_leaked_moles * leak_to_tropo_fraction
moles_mixed_to_strato_per_year = (total_leaked_moles * leak_to_strato_fraction) / years_for_strato_mix

# --- 4. Simulate the fate of the leaked methane over 3 years ---
# Initialize the amount of *additional* (leaked) methane in each layer
leaked_ch4_in_tropo = 0.0
leaked_ch4_in_strato = 0.0

# --- Year 1 ---
# Mixing
leaked_ch4_in_tropo += moles_mixed_to_tropo_y1
leaked_ch4_in_strato += moles_mixed_to_strato_per_year
# Oxidation
# 5% reduction in the troposphere, 0% in the stratosphere
leaked_ch4_in_tropo *= (1 - tropo_oxidation_rate_y1)

# --- Year 2 ---
# Mixing
leaked_ch4_in_strato += moles_mixed_to_strato_per_year
# Oxidation
# 3% reduction in both layers
leaked_ch4_in_tropo *= (1 - annual_oxidation_rate_subsequent)
leaked_ch4_in_strato *= (1 - annual_oxidation_rate_subsequent)

# --- Year 3 ---
# Mixing
leaked_ch4_in_strato += moles_mixed_to_strato_per_year
# Oxidation
# 3% reduction in both layers
leaked_ch4_in_tropo *= (1 - annual_oxidation_rate_subsequent)
leaked_ch4_in_strato *= (1 - annual_oxidation_rate_subsequent)

# Final remaining moles after 3 years
final_remaining_moles_leaked = leaked_ch4_in_tropo + leaked_ch4_in_strato

# --- 5. Calculate total moles of air in the atmosphere ---
total_mass_atmosphere_g = total_mass_atmosphere_kg * 1e3
total_moles_air = total_mass_atmosphere_g / mw_air

# --- 6. Calculate the final concentration increase in ppb ---
concentration_increase_ppb = (final_remaining_moles_leaked / total_moles_air) * 1e9

# --- 7. Print the results and the final equation ---
print("After 3 years of mixing and oxidation:")
print(f"Remaining moles from leak in Troposphere: {leaked_ch4_in_tropo:.4e} moles")
print(f"Remaining moles from leak in Stratosphere: {leaked_ch4_in_strato:.4e} moles")
print(f"Total remaining moles of leaked CH4: {final_remaining_moles_leaked:.4e} moles")
print(f"Total moles of air in atmosphere: {total_moles_air:.4e} moles")
print("\nFinal Equation for Concentration Increase (ppb):")
print(f"( {final_remaining_moles_leaked:.4e} moles_CH4 / {total_moles_air:.4e} moles_air ) * 1,000,000,000 = {concentration_increase_ppb:.4f} ppb")

# Capture the output
output = sys.stdout.getvalue()
# Restore stdout
sys.stdout = stdout_backup
# Print the captured output
print(output)
print(f"<<<{concentration_increase_ppb:.4f}>>>")