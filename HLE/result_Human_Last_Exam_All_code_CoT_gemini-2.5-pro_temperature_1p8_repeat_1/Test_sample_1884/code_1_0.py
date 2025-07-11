import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years.
    """
    # --- Given Constants ---
    ch4_release_metric_tons = 250000
    molar_mass_ch4 = 16.0  # g/mol
    total_atmosphere_mass_kg = 5.1e18  # kg
    
    # --- Standard Constants ---
    molar_mass_air_approx = 29.0 # g/mol, standard approximation for dry air
    
    # --- Calculations ---
    
    # Step 1: Calculate initial moles of CH4 released
    print("Step 1: Calculate the initial moles of methane (CH4) released.")
    ch4_release_g = ch4_release_metric_tons * 1e6
    initial_moles_ch4 = ch4_release_g / molar_mass_ch4
    print(f"Total CH4 released = ({ch4_release_metric_tons} metric tons * 1,000,000 g/ton) / {molar_mass_ch4} g/mol = {initial_moles_ch4:.4e} moles\n")

    # Step 2: Track methane amount over 3 years
    # Year 1: 80% is in the troposphere and reduced by 5%, 20% is separate
    print("Step 2: Calculate the moles of CH4 remaining after Year 1.")
    moles_in_troposphere_y1 = initial_moles_ch4 * 0.80
    moles_for_stratosphere = initial_moles_ch4 * 0.20
    reduction_troposphere_y1 = 0.05
    
    moles_after_y1 = (moles_in_troposphere_y1 * (1 - reduction_troposphere_y1)) + moles_for_stratosphere
    print(f"Year 1 Moles = (({initial_moles_ch4:.4e} * 0.80) * (1 - {reduction_troposphere_y1})) + ({initial_moles_ch4:.4e} * 0.20) = {moles_after_y1:.4e} moles\n")

    # Year 2: Total amount from year 1 is reduced by an additional 3%
    print("Step 3: Calculate the moles of CH4 remaining after Year 2.")
    reduction_y2 = 0.03
    moles_after_y2 = moles_after_y1 * (1 - reduction_y2)
    print(f"Year 2 Moles = {moles_after_y1:.4e} * (1 - {reduction_y2}) = {moles_after_y2:.4e} moles\n")

    # Year 3: Total amount from year 2 is reduced by another 3%
    print("Step 4: Calculate the moles of CH4 remaining after Year 3.")
    reduction_y3 = 0.03
    final_moles_ch4 = moles_after_y2 * (1 - reduction_y3)
    print(f"Year 3 Moles (Final) = {moles_after_y2:.4e} * (1 - {reduction_y3}) = {final_moles_ch4:.4e} moles\n")

    # Step 5: Calculate total moles of air in the atmosphere
    print("Step 5: Calculate the total moles of air in the atmosphere.")
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1e3
    total_moles_air = total_atmosphere_mass_g / molar_mass_air_approx
    print(f"Total Moles of Air = ({total_atmosphere_mass_kg:.1e} kg * 1000 g/kg) / {molar_mass_air_approx} g/mol = {total_moles_air:.4e} moles\n")

    # Step 6: Calculate final concentration increase in ppb
    print("Step 6: Calculate the final concentration increase in parts per billion (ppb).")
    ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9
    print(f"Increase (ppb) = (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
    print(f"Increase (ppb) = ({final_moles_ch4:.4e} / {total_moles_air:.4e}) * 1,000,000,000")
    print(f"Final Methane Increase = {ppb_increase:.4f} ppb")
    
    return ppb_increase

if __name__ == '__main__':
    final_answer = calculate_methane_increase()
    print(f"\n<<<The final calculated increase in atmospheric methane concentration is {final_answer:.4f} ppb.>>>")
    print(f"<<<{final_answer:.4f}>>>")
