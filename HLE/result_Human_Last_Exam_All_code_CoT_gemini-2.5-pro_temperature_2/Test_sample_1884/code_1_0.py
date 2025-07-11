import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Calculation starts here ---

# Step 1: Define constants from the problem description
initial_ch4_mass_tons = 250000.0
# Convert initial mass from metric tons to kilograms
initial_ch4_mass_kg = initial_ch4_mass_tons * 1000.0

total_atmosphere_mass_kg = 5.1 * 10**18
molar_mass_ch4_g_mol = 16.0
# The average molar mass of Earth's air is approximately 29 g/mol
molar_mass_air_g_mol = 29.0

# Define the oxidation parameters
tropospheric_mix_fraction = 0.80
y1_oxidation_rate = 0.05
subsequent_oxidation_rate = 0.03

# Step 2: Calculate the methane mass remaining after 3 years of oxidation

# In Year 1, oxidation of 5% occurs only on the 80% of methane in the troposphere.
# The total reduction is relative to the initial mass.
y1_reduction_factor = tropospheric_mix_fraction * y1_oxidation_rate
mass_after_y1 = initial_ch4_mass_kg * (1 - y1_reduction_factor)

# In Year 2, a 3% oxidation occurs on the remaining mass.
mass_after_y2 = mass_after_y1 * (1 - subsequent_oxidation_rate)

# In Year 3, another 3% oxidation occurs on the new remaining mass.
final_ch4_mass_kg = mass_after_y2 * (1 - subsequent_oxidation_rate)

# Step 3: Calculate the final concentration increase in ppb

# The ppb concentration is a mole (or volume) ratio. We convert the mass ratio to a mole ratio
# using the molar masses of air and methane.
# ppb = (mass_ratio) * (M_air / M_ch4) * 1e9
molar_mass_conversion = molar_mass_air_g_mol / molar_mass_ch4_g_mol
final_ppb_increase = (final_ch4_mass_kg / total_atmosphere_mass_kg) * molar_mass_conversion * 1e9

# Step 4: Print the process and the final result as requested.
print("Step-by-step Calculation of Final Methane Mass:")
print(f"Initial CH₄ mass = {initial_ch4_mass_kg:,.0f} kg")
print(f"Mass after Year 1 oxidation = {initial_ch4_mass_kg:,.0f} * (1 - {tropospheric_mix_fraction} * {y1_oxidation_rate}) = {mass_after_y1:,.0f} kg")
print(f"Mass after Year 2 oxidation = {mass_after_y1:,.0f} * (1 - {subsequent_oxidation_rate}) = {mass_after_y2:,.0f} kg")
print(f"Final mass after Year 3 oxidation = {mass_after_y2:,.0f} * (1 - {subsequent_oxidation_rate}) = {final_ch4_mass_kg:,.0f} kg")

print("\nFinal Equation for Concentration Increase (ppb):")
print("ppb = (Final CH₄ Mass / Total Atmosphere Mass) * (Molar Mass of Air / Molar Mass of CH₄) * 10^9")
print(f"ppb = ({final_ch4_mass_kg:,.0f} kg / {total_atmosphere_mass_kg:.1e} kg) * ({molar_mass_air_g_mol} / {molar_mass_ch4_g_mol}) * 1,000,000,000")
print(f"ppb = {final_ppb_increase:.5f}")

# --- End of Calculation ---

# Restore stdout and get the output
sys.stdout = old_stdout
output_str = captured_output.getvalue()

# Print the captured output
print(output_str)

# Finally, print the answer in the specified format
print(f"<<<{final_ppb_increase:.5f}>>>")
