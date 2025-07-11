import sys

# Step 1: Define constants and initial methane mass
total_ch4_leak_tons = 250000.0
tons_to_g = 1e6
total_ch4_leak_g = total_ch4_leak_tons * tons_to_g

# Atmospheric constants
total_atmosphere_mass_kg = 5.1e18
total_atmosphere_mass_g = total_atmosphere_mass_kg * 1000.0

# Molecular weights for conversion from mass ratio to volume ratio
mw_ch4_g_per_mol = 16.0
mw_air_g_per_mol = 29.0 # Average molecular weight of air

# Separate initial mass into two portions
mass_portion_A_initial = total_ch4_leak_g * 0.80 # Enters troposphere
mass_portion_B_source = total_ch4_leak_g * 0.20 # Enters stratosphere

# Step 2: Track decay of the tropospheric portion (Portion A)
# Year 1: 5% decay
mass_A_y1 = mass_portion_A_initial * (1 - 0.05)
# Year 2: 3% decay
mass_A_y2 = mass_A_y1 * (1 - 0.03)
# Year 3: 3% decay
mass_A_final = mass_A_y2 * (1 - 0.03)

# Step 3: Track mixing and decay of the stratospheric portion (Portion B)
# It is added in three equal chunks over 3 years.
strato_chunk = mass_portion_B_source / 3.0

# Mass from Portion B at the end of each year
mass_B_y1 = strato_chunk  # Added in Y1, no decay specified for it in Y1
mass_B_y2 = (mass_B_y1 * (1 - 0.03)) + strato_chunk # Mass from Y1 decays, new chunk added
mass_B_final = (mass_B_y2 * (1 - 0.03)) + strato_chunk # Mass from Y2 decays, new chunk added

# Step 4: Calculate total remaining methane mass
total_remaining_ch4_g = mass_A_final + mass_B_final

# Step 5: Convert total remaining mass to concentration increase in ppb
# ppb_increase = (total_remaining_ch4_g / total_atmosphere_mass_g) * (mw_air_g_per_mol / mw_ch4_g_per_mol) * 1e9

# To fulfill the request "output each number in the final equation", we print the equation with its values.
# Final equation: (Total CH4 Remaining / Total Atmosphere Mass) * (MW_Air / MW_CH4) * 1e9 = ppb increase
final_calc_str = (
    f"({total_remaining_ch4_g:.3e} g / {total_atmosphere_mass_g:.3e} g) * "
    f"({mw_air_g_per_mol} / {mw_ch4_g_per_mol}) * 1e9 = "
    f"{(total_remaining_ch4_g / total_atmosphere_mass_g) * (mw_air_g_per_mol / mw_ch4_g_per_mol) * 1e9:.4f} ppb"
)
print(final_calc_str)

final_ppb_increase = (total_remaining_ch4_g / total_atmosphere_mass_g) * (mw_air_g_per_mol / mw_ch4_g_per_mol) * 1e9
sys.stdout.write(f"\n<<<{final_ppb_increase:.4f}>>>")