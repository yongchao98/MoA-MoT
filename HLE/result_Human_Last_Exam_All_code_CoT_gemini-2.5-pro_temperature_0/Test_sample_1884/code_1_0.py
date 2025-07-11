import sys

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after a 250,000 metric ton release, considering mixing and oxidation over 3 years.
    """
    # --- Given Constants ---
    total_ch4_release_tons = 250000
    molecular_weight_ch4 = 16.0  # g/mol
    total_atmosphere_mass_kg = 5.1e18
    
    # Standard value for average molar mass of dry air, needed for ppb calculation
    avg_molar_mass_air = 28.97 # g/mol

    # --- Step 1: Calculate the total mass of methane released in grams ---
    total_ch4_release_g = total_ch4_release_tons * 1e6
    print(f"Step 1: Calculate the total mass of methane released in grams.")
    print(f"Equation: {total_ch4_release_tons} metric tons * 1,000,000 g/ton = {total_ch4_release_g:.2e} g\n")

    # --- Step 2: Model the mixing and oxidation for Year 1 ---
    ch4_in_troposphere_g = total_ch4_release_g * 0.80
    ch4_unmixed_g = total_ch4_release_g * 0.20
    tropospheric_oxidation_rate = 0.05
    ch4_after_y1_oxidation = ch4_in_troposphere_g * (1 - tropospheric_oxidation_rate)
    ch4_mass_after_y1 = ch4_after_y1_oxidation + ch4_unmixed_g
    print(f"Step 2: Calculate the methane mass remaining after Year 1.")
    print(f"Equation: ({total_ch4_release_g:.2e} g * 0.80 * (1 - {tropospheric_oxidation_rate})) + ({total_ch4_release_g:.2e} g * 0.20) = {ch4_mass_after_y1:.2e} g\n")

    # --- Step 3: Calculate the methane mass remaining after Year 2 ---
    subsequent_oxidation_rate = 0.03
    ch4_mass_after_y2 = ch4_mass_after_y1 * (1 - subsequent_oxidation_rate)
    print(f"Step 3: Calculate the methane mass remaining after Year 2.")
    print(f"Equation: {ch4_mass_after_y1:.2e} g * (1 - {subsequent_oxidation_rate}) = {ch4_mass_after_y2:.2e} g\n")

    # --- Step 4: Calculate the final methane mass remaining after Year 3 ---
    ch4_mass_after_y3 = ch4_mass_after_y2 * (1 - subsequent_oxidation_rate)
    print(f"Step 4: Calculate the final methane mass remaining after Year 3.")
    print(f"Equation: {ch4_mass_after_y2:.2e} g * (1 - {subsequent_oxidation_rate}) = {ch4_mass_after_y3:.2e} g\n")

    # --- Step 5: Convert the final methane mass to moles ---
    final_ch4_moles = ch4_mass_after_y3 / molecular_weight_ch4
    print(f"Step 5: Convert the final methane mass to moles.")
    print(f"Equation: {ch4_mass_after_y3:.2e} g / {molecular_weight_ch4} g/mol = {final_ch4_moles:.2e} moles\n")

    # --- Step 6: Calculate the total moles of air in the atmosphere ---
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1000
    total_air_moles = total_atmosphere_mass_g / avg_molar_mass_air
    print(f"Step 6: Calculate the total moles of air in the atmosphere.")
    print(f"Equation: ({total_atmosphere_mass_kg:.1e} kg * 1000 g/kg) / {avg_molar_mass_air} g/mol = {total_air_moles:.2e} moles\n")

    # --- Step 7: Calculate the final concentration increase in parts per billion (ppb) ---
    ppb_increase = (final_ch4_moles / total_air_moles) * 1e9
    print(f"Step 7: Calculate the final concentration increase in parts per billion (ppb).")
    print(f"Final Equation: ({final_ch4_moles:.2e} moles CH₄ / {total_air_moles:.2e} moles air) * 1,000,000,000 = {ppb_increase:.2f} ppb\n")
    
    # --- Final Answer ---
    # Use sys.stdout to ensure the final answer is the last thing printed
    sys.stdout.write(f"<<<{ppb_increase:.2f}>>>\n")

calculate_methane_increase()