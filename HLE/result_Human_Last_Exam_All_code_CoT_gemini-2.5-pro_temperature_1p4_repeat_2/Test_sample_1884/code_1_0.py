import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Create a new output stream
new_stdout = io.StringIO()
# Redirect stdout
sys.stdout = new_stdout

# --- Start of the student's solution ---

# Given values
total_methane_release_tons = 250000.0
troposphere_mixing_fraction = 0.8
stratosphere_mixing_fraction = 0.2
year1_tropospheric_oxidation_reduction = 0.05
subsequent_years_oxidation_reduction = 0.03
years_of_subsequent_oxidation = 2  # Year 2 and Year 3

molecular_weight_ch4 = 16.0  # g/mol
average_molecular_weight_air = 29.0  # g/mol, average molecular weight of dry air
total_mass_atmosphere_kg = 5.1 * 10**18

# Convert total release to kg
total_methane_release_kg = total_methane_release_tons * 1000

# Calculate the initial mass distributed to the troposphere and stratosphere
initial_mass_troposphere_kg = total_methane_release_kg * troposphere_mixing_fraction
initial_mass_stratosphere_kg = total_methane_release_kg * stratosphere_mixing_fraction

# Calculate the remaining mass in the troposphere after 3 years of oxidation
# Year 1: 5% reduction
# Years 2 & 3: 3% reduction each year
remaining_mass_troposphere_kg = initial_mass_troposphere_kg * (1 - year1_tropospheric_oxidation_reduction) * ((1 - subsequent_years_oxidation_reduction) ** years_of_subsequent_oxidation)

# Calculate the remaining mass in the stratosphere after 3 years of oxidation
# No reduction in Year 1
# Years 2 & 3: 3% reduction each year
remaining_mass_stratosphere_kg = initial_mass_stratosphere_kg * ((1 - subsequent_years_oxidation_reduction) ** years_of_subsequent_oxidation)

# Calculate the total remaining mass of methane
total_remaining_methane_kg = remaining_mass_troposphere_kg + remaining_mass_stratosphere_kg

# Calculate the final mass fraction in the total atmosphere
final_mass_fraction = total_remaining_methane_kg / total_mass_atmosphere_kg

# Convert mass fraction to mole fraction (ppb is a mole-based measure)
# Mole Fraction = Mass Fraction * (MW_air / MW_substance)
final_mole_fraction = final_mass_fraction * (average_molecular_weight_air / molecular_weight_ch4)

# Convert mole fraction to parts per billion (ppb)
final_ppb_increase = final_mole_fraction * 1e9

# Print the final equation with all the numbers
print("Equation to calculate the final ppb increase:")
print(f"Final ppb = ((({initial_mass_troposphere_kg:.1e} kg * (1 - {year1_tropospheric_oxidation_reduction}) * (1 - {subsequent_years_oxidation_reduction})^2) + "
      f"({initial_mass_stratosphere_kg:.1e} kg * (1 - {subsequent_years_oxidation_reduction})^2)) / "
      f"{total_mass_atmosphere_kg:.1e} kg) * ({average_molecular_weight_air} / {molecular_weight_ch4}) * 1,000,000,000")

# Print the final result
print(f"\nFinal increase in atmospheric methane concentration: {final_ppb_increase:.5f} ppb")

# --- End of the student's solution ---

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the new_stdout
output = new_stdout.getvalue()

# Print the captured output
print(output)
print(f"<<<{final_ppb_increase:.5f}>>>")