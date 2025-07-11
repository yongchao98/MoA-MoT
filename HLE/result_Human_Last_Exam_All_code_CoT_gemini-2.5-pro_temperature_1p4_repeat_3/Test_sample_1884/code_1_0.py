import sys

# Plan:
# 1. Calculate the mass of methane remaining after 3 years of oxidation.
# 2. Calculate the mass mixing ratio in the total atmosphere.
# 3. Convert the mass mixing ratio to a molar mixing ratio (ppb).

# --- Given Constants ---
initial_methane_release_tons = 250000  # metric tons
mass_atmosphere_kg = 5.1 * 10**18     # kg
mw_ch4 = 16.0                         # g/mol
mw_air = 29.0                         # g/mol (average molecular weight of air)
year1_oxidation_rate = 0.05           # 5%
subsequent_oxidation_rate = 0.03      # 3%

# --- Step 1: Calculate Methane Decay ---

# Mass after 1 year (5% reduction)
mass_after_1_year = initial_methane_release_tons * (1 - year1_oxidation_rate)

# Mass after 2 years (additional 3% reduction)
mass_after_2_years = mass_after_1_year * (1 - subsequent_oxidation_rate)

# Mass after 3 years (final 3% reduction)
final_mass_ch4_tons = mass_after_2_years * (1 - subsequent_oxidation_rate)

# --- Step 2: Calculate Final Mass Mixing Ratio ---

# Convert final mass from metric tons to kilograms
final_mass_ch4_kg = final_mass_ch4_tons * 1000

# Mass mixing ratio (kg CH4 / kg air)
mass_mixing_ratio = final_mass_ch4_kg / mass_atmosphere_kg

# --- Step 3 & 4: Convert to Molar Ratio (ppb) ---

# Molar mixing ratio = mass_mixing_ratio * (mw_air / mw_ch4)
# ppb = molar_mixing_ratio * 1e9
ppb_increase = mass_mixing_ratio * (mw_air / mw_ch4) * 1e9

# --- Output the result ---
print("This script calculates the final increase in atmospheric methane concentration in ppb after 3 years.")
print("-" * 50)
print(f"Initial Methane Release: {initial_methane_release_tons} metric tons")
print(f"Methane remaining after 3 years: {final_mass_ch4_tons:.2f} metric tons")
print(f"Final methane mass in kg: {final_mass_ch4_kg:.2e} kg")
print("-" * 50)
print("Final Calculation:")
print("ppb increase = (Final CH4 Mass / Total Atmosphere Mass) * (MW_air / MW_CH4) * 1e9")
print(f"ppb increase = ({final_mass_ch4_kg:.4e} kg / {mass_atmosphere_kg:.1e} kg) * ({mw_air} / {mw_ch4}) * 1e9")
print(f"Result: The increase in atmospheric methane concentration after 3 years is {ppb_increase:.2f} ppb.")

# The final result in the specified format
sys.stdout.write(f"\n<<<{ppb_increase:.2f}>>>")