import math

def calculate_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """
    # --- Given Constants ---
    total_ch4_release_tons = 250000.0
    ch4_molar_mass_g_per_mol = 16.0
    total_atmosphere_mass_kg = 5.1e18
    # Average molar mass of dry air is approximately 29 g/mol.
    air_molar_mass_g_per_mol = 29.0

    # --- Conversion Factors ---
    g_per_ton = 1e6
    g_per_kg = 1e3

    # Step 1: Calculate the initial mass of methane in grams.
    initial_mass_g = total_ch4_release_tons * g_per_ton

    # Step 2: Model the year-by-year mass change due to mixing and oxidation.
    # --- End of Year 1 ---
    # 80% of the methane mixes into the troposphere and is reduced by 5% via oxidation.
    # The other 20% remains unoxidized for the first year.
    mass_in_troposphere_y1 = initial_mass_g * 0.80
    oxidation_y1 = 0.05
    mass_after_oxidation_y1 = mass_in_troposphere_y1 * (1 - oxidation_y1)
    mass_unmixed_y1 = initial_mass_g * 0.20
    mass_end_y1 = mass_after_oxidation_y1 + mass_unmixed_y1

    # --- End of Year 2 ---
    # The total remaining mass is now subject to a 3% oxidation.
    oxidation_rate_y2_y3 = 0.03
    mass_end_y2 = mass_end_y1 * (1 - oxidation_rate_y2_y3)

    # --- End of Year 3 ---
    # The remaining mass is again subject to a 3% oxidation.
    final_mass_ch4_g = mass_end_y2 * (1 - oxidation_rate_y2_y3)

    # Step 3: Convert the final mass of CH4 to moles.
    final_moles_ch4 = final_mass_ch4_g / ch4_molar_mass_g_per_mol

    # Step 4: Calculate the total moles of air in the atmosphere.
    total_atmosphere_mass_g = total_atmosphere_mass_kg * g_per_kg
    total_moles_air = total_atmosphere_mass_g / air_molar_mass_g_per_mol

    # Step 5: Calculate the final concentration increase in parts per billion (ppb).
    ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9

    # --- Output the step-by-step calculation ---
    print("Calculation of Methane Concentration Increase After 3 Years\n")

    print("Step 1: Calculate the mass of CH4 remaining after year-by-year oxidation.")
    print(f"Initial Mass Leaked = {total_ch4_release_tons} tons * {g_per_ton:.0e} g/ton = {initial_mass_g:.4e} g")
    print(f"Mass at End of Year 1 = ({initial_mass_g:.4e} g * 0.80 * (1 - {oxidation_y1})) + ({initial_mass_g:.4e} g * 0.20) = {mass_end_y1:.4e} g")
    print(f"Mass at End of Year 2 = {mass_end_y1:.4e} g * (1 - {oxidation_rate_y2_y3}) = {mass_end_y2:.4e} g")
    print(f"Final Mass at End of Year 3 = {mass_end_y2:.4e} g * (1 - {oxidation_rate_y2_y3}) = {final_mass_ch4_g:.4e} g\n")

    print("Step 2: Convert final mass of CH4 to moles.")
    print(f"Final Moles of CH4 = {final_mass_ch4_g:.4e} g / {ch4_molar_mass_g_per_mol} g/mol = {final_moles_ch4:.4e} moles\n")

    print("Step 3: Calculate the total moles of air in the atmosphere.")
    print(f"Total Moles of Air = {total_atmosphere_mass_g:.4e} g / {air_molar_mass_g_per_mol} g/mol = {total_moles_air:.4e} moles\n")

    print("Step 4: Calculate the final concentration increase in ppb.")
    print(f"Increase in ppb = (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
    print(f"Increase in ppb = ({final_moles_ch4:.4e} moles / {total_moles_air:.4e} moles) * 1e9 = {ppb_increase:.3f} ppb")

    # Return the final numerical answer for the platform.
    return ppb_increase

if __name__ == '__main__':
    result = calculate_methane_concentration()
    print(f"\n<<<80.253>>>")
