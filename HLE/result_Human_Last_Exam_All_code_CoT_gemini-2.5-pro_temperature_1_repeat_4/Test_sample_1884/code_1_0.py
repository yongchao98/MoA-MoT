import sys

# --- Constants from the problem statement ---
# Initial mass of methane (CH4) released in metric tons
initial_ch4_mass_tons = 250000
# Molar mass of CH4 in g/mol
molar_mass_ch4_g_mol = 16.0
# Total mass of the atmosphere in kg
total_atmo_mass_kg = 5.1e18
# Average molar mass of air in g/mol (approx. 29 g/mol is a standard value)
molar_mass_air_g_mol = 29.0

# --- Proportions and Rates ---
# Fraction of CH4 that mixes into the troposphere in year 1
frac_tropo = 0.8
# Oxidation rate in the troposphere in year 1
ox_rate_y1 = 0.05
# Oxidation rate in the whole atmosphere for subsequent years
ox_rate_y23 = 0.03

# --- Step 1: Convert initial mass to kg ---
initial_ch4_mass_kg = initial_ch4_mass_tons * 1000

# --- Step 2: Calculate remaining mass after Year 1 ---
# Mass that enters the troposphere
mass_in_tropo_y1 = initial_ch4_mass_kg * frac_tropo
# Mass remaining from the tropospheric portion after 5% oxidation
mass_from_tropo_y1_end = mass_in_tropo_y1 * (1 - ox_rate_y1)
# Mass that was held back and mixes in later
mass_held_back = initial_ch4_mass_kg * (1 - frac_tropo)

# --- Step 3: Calculate remaining mass after Years 2 and 3 ---
# Total mass at the start of Year 2 (fully mixed)
total_mass_y2_start = mass_from_tropo_y1_end + mass_held_back
# Mass remaining after 3% oxidation in Year 2
total_mass_y2_end = total_mass_y2_start * (1 - ox_rate_y23)
# Final mass remaining after 3% oxidation in Year 3
final_ch4_mass_kg = total_mass_y2_end * (1 - ox_rate_y23)

# --- Step 4: Calculate the final concentration increase in ppb ---
# Convert final CH4 mass to moles
final_ch4_moles = (final_ch4_mass_kg * 1000) / molar_mass_ch4_g_mol
# Calculate total moles of air in the atmosphere
total_air_moles = (total_atmo_mass_kg * 1000) / molar_mass_air_g_mol
# Calculate the final increase in ppb
ppb_increase = (final_ch4_moles / total_air_moles) * 1e9

# --- Final Output ---
# Print the complete equation with all the numbers, as requested.
# Note: The expression is evaluated to show the final result matches the step-by-step calculation.
# It represents: ( ( (initial_mass * frac_tropo * (1-ox1)) + (initial_mass * (1-frac_tropo)) ) * (1-ox23) * (1-ox23) * 1000 / molar_mass_ch4 ) / ( total_atmo_mass * 1000 / molar_mass_air ) * 1e9
print("Final Equation for Methane Increase (in ppb):")
print(
    f"( ( ( {initial_ch4_mass_kg} * {frac_tropo} * {1-ox_rate_y1} ) + ( {initial_ch4_mass_kg} * {1-frac_tropo} ) ) * {1-ox_rate_y23} * {1-ox_rate_y23} * 1000 / {molar_mass_ch4_g_mol} ) / "
    f"( {total_atmo_mass_kg} * 1000 / {molar_mass_air_g_mol} ) * 1000000000 = {ppb_increase:.5f}"
)

# Suppress the thought process from being part of the final answer block
if 'ipykernel' in sys.modules:
    pass # In a notebook environment, just let the print happen.
else:
    # In other environments, this will be the only output besides the print statements.
    sys.stdout = open(os.devnull, 'w')
    
# The final answer in the required format.
final_answer = f"<<<{ppb_increase:.5f}>>>"