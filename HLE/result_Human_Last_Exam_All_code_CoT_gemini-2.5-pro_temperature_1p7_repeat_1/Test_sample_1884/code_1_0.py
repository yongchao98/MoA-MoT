import math

def calculate_methane_increase():
    """
    This function calculates the increase in atmospheric methane concentration in ppb
    after a 250,000 metric ton release, considering mixing and oxidation over 3 years.
    """

    # --- Step 1: Define Constants from the problem statement ---
    initial_ch4_release_tons = 250000.0  # metric tons
    total_atmosphere_mass_kg = 5.1 * 10**18  # kg
    molar_mass_ch4_g_mol = 16.0  # g/mol
    molar_mass_air_g_mol = 29.0  # g/mol (average for air)
    
    # Mixing and Oxidation Parameters
    tropospheric_fraction = 0.80
    stratospheric_fraction = 0.20
    oxidation_loss_y1_rate = 0.05  # 5%
    oxidation_loss_subsequent_rate = 0.03 # 3%

    # --- Step 2: Calculate the fate of methane mass over 3 years ---
    initial_ch4_release_kg = initial_ch4_release_tons * 1000

    print("--- Methane Mass Calculation Over 3 Years ---")
    print(f"Initial methane release: {initial_ch4_release_kg:,.0f} kg")

    # End of Year 1: 80% is oxidized by 5%, 20% is untouched.
    mass_end_y1 = initial_ch4_release_kg * (tropospheric_fraction * (1 - oxidation_loss_y1_rate) + stratospheric_fraction)
    print(f"Mass remaining after Year 1: {mass_end_y1:,.0f} kg")

    # End of Year 2: The entire remaining mass is reduced by 3%.
    mass_end_y2 = mass_end_y1 * (1 - oxidation_loss_subsequent_rate)
    print(f"Mass remaining after Year 2: {mass_end_y2:,.0f} kg")

    # End of Year 3: The remaining mass is again reduced by 3%.
    final_ch4_mass_increase_kg = mass_end_y2 * (1 - oxidation_loss_subsequent_rate)
    print(f"Final methane mass remaining after 3 years: {final_ch4_mass_increase_kg:,.0f} kg\n")

    # --- Step 3: Convert the final mass to a concentration increase in ppb ---
    # The conversion formula from a mass ratio to a mole ratio (ppb) is:
    # ppb = (Mass CH4 / Mass Atmosphere) * (Molar Mass of Air / Molar Mass of CH4) * 10^9
    ppb_increase = (final_ch4_mass_increase_kg / total_atmosphere_mass_kg) * (molar_mass_air_g_mol / molar_mass_ch4_g_mol) * 1e9

    print("--- Final Concentration Calculation ---")
    print("The final concentration increase is calculated by converting the final mass to a mole-based concentration (ppb).\n")
    
    print("Final Equation:")
    print(f"ppb_increase = ({final_ch4_mass_increase_kg:,.2f} kg / {total_atmosphere_mass_kg:.1e} kg) * ({molar_mass_air_g_mol} g/mol / {molar_mass_ch4_g_mol} g/mol) * 1,000,000,000")
    
    return ppb_increase

# Execute the calculation and print the final answer in the required format.
final_answer = calculate_methane_increase()
print(f"<<<{final_answer:.2f}>>>")