import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# Define constants based on the problem statement
total_release_tons = 250000  # metric tons
atmosphere_mass_kg = 5.1e18
molar_mass_ch4 = 16.0  # g/mol
# Use a standard value for the average molar mass of dry air
molar_mass_air = 28.97  # g/mol

# --- Step 1: Calculate initial mass in kg ---
total_release_kg = total_release_tons * 1000

# --- Step 2: Track mass reduction over 3 years ---
# At the start of Year 1
mass_in_tropo_y1 = total_release_kg * 0.80
mass_to_be_mixed = total_release_kg * 0.20

# At the end of Year 1, after 5% oxidation in the troposphere
mass_in_tropo_end_y1 = mass_in_tropo_y1 * (1 - 0.05)
total_mass_end_y1 = mass_in_tropo_end_y1 + mass_to_be_mixed

# At the end of Year 2, after a 3% reduction on the total remaining amount
total_mass_end_y2 = total_mass_end_y1 * (1 - 0.03)

# At the end of Year 3, after another 3% reduction
final_ch4_mass_kg = total_mass_end_y2 * (1 - 0.03)

# --- Step 3: Calculate final concentration increase in ppb ---
# ppb (by mole) = (mass_fraction) * (molar_mass_air / molar_mass_ch4) * 1e9
mass_fraction = final_ch4_mass_kg / atmosphere_mass_kg
ppb_increase = mass_fraction * (molar_mass_air / molar_mass_ch4) * 1e9

# --- Step 4: Print the detailed calculation ---
print("Step 1: Calculate the initial mass of methane released.")
print(f"Total methane released = {total_release_tons} metric tons = {total_release_kg:.2e} kg\n")

print("Step 2: Calculate the mass of methane remaining after 3 years.")
print(f"Mass at the end of Year 1 after 5% oxidation in the troposphere:")
print(f"  (({total_release_kg:.2e} kg * 0.80) * (1 - 0.05)) + ({total_release_kg:.2e} kg * 0.20) = {total_mass_end_y1:.3e} kg")
print(f"Mass at the end of Year 2 after a further 3% oxidation:")
print(f"  {total_mass_end_y1:.3e} kg * (1 - 0.03) = {total_mass_end_y2:.3e} kg")
print(f"Final mass at the end of Year 3 after a final 3% oxidation:")
print(f"  {total_mass_end_y2:.3e} kg * (1 - 0.03) = {final_ch4_mass_kg:.3e} kg\n")

print("Step 3: Calculate the final concentration increase in ppb.")
print("The concentration in ppb is calculated as: (Mass CH4 / Mass Atmosphere) * (Molar Mass Air / Molar Mass CH4) * 10^9\n")

print("Final Equation:")
print(f"Increase (ppb) = ({final_ch4_mass_kg:.5e} / {atmosphere_mass_kg:.5e}) * ({molar_mass_air} / {molar_mass_ch4}) * 10^9")

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the buffer
output = buffer.getvalue()
print(output)

print(f"The calculated increase in atmospheric methane concentration after 3 years is {ppb_increase:.5f} ppb.")

print(f"\n<<<{ppb_increase:.5f}>>>")