import sys

# Step 1: Define constants and initial values
ch4_emission_tons = 250000.0  # metric tons
total_atmosphere_mass_kg = 5.1e18  # kg
mw_ch4 = 16.0  # g/mol
mw_air = 29.0  # Average molecular weight of air, g/mol

# Unit conversions for consistency
ch4_emission_g = ch4_emission_tons * 1e6  # convert tons to grams
total_atmosphere_mass_g = total_atmosphere_mass_kg * 1e3  # convert kg to grams

# Methane distribution
ch4_in_troposphere_initial = ch4_emission_g * 0.80
ch4_for_stratosphere = ch4_emission_g * 0.20

# Step 2: Simulate Year 1
# 5% oxidation of the methane in the troposphere
ch4_after_year1 = ch4_in_troposphere_initial * (1 - 0.05)

# Step 3: Simulate Year 2
# The remaining 20% mixes into the atmosphere
total_ch4_start_year2 = ch4_after_year1 + ch4_for_stratosphere
# 3% oxidation of the total amount
ch4_after_year2 = total_ch4_start_year2 * (1 - 0.03)

# Step 4: Simulate Year 3
# 3% oxidation of the remaining amount
final_ch4_mass_g = ch4_after_year2 * (1 - 0.03)

# Step 5: Calculate the final concentration increase in ppb
# ppb = (mass_fraction) * (molar_mass_ratio) * 1e9
ppb_increase = (final_ch4_mass_g / total_atmosphere_mass_g) * (mw_air / mw_ch4) * 1e9

# Step 6: Output the final equation and result
print("Calculation for the final concentration increase:")
# Use f-strings to format the output string with the calculated values
# The format specifier .4e prints the number in scientific notation with 4 decimal places
print(
    f"ppb = ({final_ch4_mass_g:.4e} g CH₄ / {total_atmosphere_mass_g:.4e} g Air) * "
    f"({mw_air} g/mol Air / {mw_ch4} g/mol CH₄) * 1e9"
)
print(f"\nFinal increase in methane concentration after 3 years: {ppb_increase:.5f} ppb")

# Required final answer format
# Using sys.stdout.write to avoid adding a newline, as per some platform requirements
sys.stdout.write(f"<<<{ppb_increase:.5f}>>>")