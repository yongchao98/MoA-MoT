import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of calculations ---

# Step 1: Define constants from the problem statement
total_ch4_release_tons = 250000.0  # metric tons
ch4_molecular_weight = 16.0  # g/mol
atmosphere_total_mass_kg = 5.1e18  # kg
air_avg_molecular_weight = 29.0  # g/mol (standard approximation for air)

troposphere_mix_fraction = 0.80
year1_oxidation_rate = 0.05
subsequent_year_oxidation_rate = 0.03

# Step 2: Calculate the initial total mass of methane in grams
total_ch4_mass_g = total_ch4_release_tons * 1e6 # 1e6 grams per metric ton

# Step 3: Partition the mass and track its reduction over 3 years
# Initial partition
mass_in_troposphere_pool = total_ch4_mass_g * troposphere_mix_fraction
mass_in_stratosphere_pool = total_ch4_mass_g * (1 - troposphere_mix_fraction)

# End of Year 1: 5% oxidation on the tropospheric pool only
mass_troposphere_pool_y1_end = mass_in_troposphere_pool * (1 - year1_oxidation_rate)
total_mass_y1_end = mass_troposphere_pool_y1_end + mass_in_stratosphere_pool

# End of Year 2: 3% oxidation on the total remaining mass
total_mass_y2_end = total_mass_y1_end * (1 - subsequent_year_oxidation_rate)

# End of Year 3: 3% oxidation on the total remaining mass
final_mass_y3_end = total_mass_y2_end * (1 - subsequent_year_oxidation_rate)

# Step 4: Convert final methane mass to moles
final_ch4_moles = final_mass_y3_end / ch4_molecular_weight

# Step 5: Calculate total moles of air in the atmosphere
atmosphere_total_mass_g = atmosphere_total_mass_kg * 1000
total_air_moles = atmosphere_total_mass_g / air_avg_molecular_weight

# Step 6: Calculate the final concentration increase in parts per billion (ppb)
final_concentration_ppb = (final_ch4_moles / total_air_moles) * 1e9

# --- Final Output Generation ---
# The prompt requires printing each number in the final equation.

print("Final Calculation Steps:")
print(f"The final mass of CH4 after 3 years is {final_mass_y3_end:.2f} g.")
print(f"Converting to moles: {final_mass_y3_end:.2f} g / {ch4_molecular_weight} g/mol = {final_ch4_moles:.2f} moles of CH4.")
print(f"The total moles of air in the atmosphere is {atmosphere_total_mass_g:.2e} g / {air_avg_molecular_weight} g/mol = {total_air_moles:.2e} moles.")
print("-" * 30)
print("Final Concentration Equation:")
print(f"Increase in CH4 (ppb) = (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
print(f"Increase in CH4 (ppb) = ({final_ch4_moles:.2f} / {total_air_moles:.2e}) * 1,000,000,000")
print(f"Result: {final_concentration_ppb:.4f} ppb")

# --- End of calculations ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)
# Finally, provide the answer in the specified format
print(f"<<<{final_concentration_ppb:.4f}>>>")