import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Create a new StringIO object
# This is a trick to get the final answer without the user having to copy it.
string_io = io.StringIO()
# Redirect stdout to the StringIO object
sys.stdout = string_io

# --- Define Constants ---
# Total methane release in metric tons and grams
total_ch4_release_tons = 250000
total_ch4_release_g = float(total_ch4_release_tons * 1e6)

# Mass of the total atmosphere in kg and grams
mass_atmosphere_kg = 5.1e18
mass_atmosphere_g = float(mass_atmosphere_kg * 1e3)

# Molar masses
molar_mass_ch4 = 16.0  # g/mol
molar_mass_air = 29.0   # Average molar mass of air in g/mol

# Oxidation rates
oxidation_y1_tropo = 0.05
oxidation_y2_y3 = 0.03

# --- Step 1: Initial Mass Allocation ---
ch4_to_troposphere = total_ch4_release_g * 0.80
ch4_to_stratosphere_total = total_ch4_release_g * 0.20
ch4_stratosphere_annual_addition = ch4_to_stratosphere_total / 3.0

print("--- Calculation Steps ---")
print(f"1. Methane Mass Allocation:")
print(f"   - Mass for Troposphere (80%): {ch4_to_troposphere:.4e} g")
print(f"   - Annual Mass for Stratosphere (20% over 3 years): {ch4_stratosphere_annual_addition:.4e} g\n")


# --- Step 2: Year-by-Year Simulation ---
mass_in_troposphere = 0.0
mass_in_stratosphere = 0.0

# Year 1
print("2. Year 1 Simulation:")
mass_in_troposphere += ch4_to_troposphere
mass_in_stratosphere += ch4_stratosphere_annual_addition
print(f"   - Added to layers: Troposphere={mass_in_troposphere:.4e} g, Stratosphere={mass_in_stratosphere:.4e} g")
mass_in_troposphere *= (1 - oxidation_y1_tropo)
print(f"   - After 5% tropospheric oxidation: Troposphere={mass_in_troposphere:.4e} g\n")

# Year 2
print("3. Year 2 Simulation:")
mass_in_stratosphere += ch4_stratosphere_annual_addition
print(f"   - After adding to stratosphere: Stratosphere={mass_in_stratosphere:.4e} g")
mass_in_troposphere *= (1 - oxidation_y2_y3)
mass_in_stratosphere *= (1 - oxidation_y2_y3)
print(f"   - After 3% oxidation in both layers: Troposphere={mass_in_troposphere:.4e} g, Stratosphere={mass_in_stratosphere:.4e} g\n")

# Year 3
print("4. Year 3 Simulation:")
mass_in_stratosphere += ch4_stratosphere_annual_addition
print(f"   - After adding to stratosphere: Stratosphere={mass_in_stratosphere:.4e} g")
mass_in_troposphere *= (1 - oxidation_y2_y3)
mass_in_stratosphere *= (1 - oxidation_y2_y3)
print(f"   - After 3% oxidation in both layers: Troposphere={mass_in_troposphere:.4e} g, Stratosphere={mass_in_stratosphere:.4e} g\n")


# --- Step 3: Final Calculation ---
total_remaining_ch4_mass = mass_in_troposphere + mass_in_stratosphere

# Convert final mass fraction to ppb by volume
final_ppb_increase = (total_remaining_ch4_mass / mass_atmosphere_g) * (molar_mass_air / molar_mass_ch4) * 1e9

print("5. Final Concentration Calculation:")
print(f"   - Total remaining CH4 mass = {mass_in_troposphere:.4e} g + {mass_in_stratosphere:.4e} g = {total_remaining_ch4_mass:.4e} g")
print("\n--- Final Equation ---")
# The final equation requires all numbers to be outputted.
print(f"Increase in CH4 (ppb) = (Total Remaining Mass / Total Atmosphere Mass) * (Molar Mass Air / Molar Mass CH4) * 1,000,000,000")
print(f"Increase in CH4 (ppb) = ({total_remaining_ch4_mass:.4e} g / {mass_atmosphere_g:.2e} g) * ({molar_mass_air} g/mol / {molar_mass_ch4} g/mol) * 1e9")
print(f"Increase in CH4 (ppb) = {final_ppb_increase:.5f}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the content from the StringIO object
output = string_io.getvalue()
# Print the captured output to the real console
print(output)
# Print the final answer in the desired format.
print(f'<<<{final_ppb_increase:.5f}>>>')