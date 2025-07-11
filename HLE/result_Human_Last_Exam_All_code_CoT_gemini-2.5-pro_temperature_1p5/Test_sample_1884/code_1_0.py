import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years following the Nord Stream leak, based on the provided data.
    """
    # Step 1: Define constants and initial values from the problem description
    methane_release_tons = 250000.0  # metric tons
    molecular_weight_ch4 = 16.0  # g/mol
    total_atmosphere_mass_kg = 5.1e18  # kg
    avg_molar_mass_air = 29.0  # g/mol, a standard approximation

    # Oxidation and distribution parameters
    troposphere_fraction = 0.80
    year1_oxidation_rate = 0.05
    subsequent_years_oxidation_rate = 0.03

    # Step 2: Calculate the total moles of methane released
    methane_release_g = methane_release_tons * 1e6  # Convert tons to grams
    total_moles_methane_released = methane_release_g / molecular_weight_ch4

    # Step 3: Model the methane decay over 3 years
    # Initial split of methane
    initial_moles_in_troposphere = total_moles_methane_released * troposphere_fraction
    initial_moles_in_stratosphere = total_moles_methane_released * (1 - troposphere_fraction)

    # After Year 1: 5% oxidation of the tropospheric part
    moles_after_y1_oxidation = initial_moles_in_troposphere * (1 - year1_oxidation_rate)
    total_moles_at_end_of_y1 = moles_after_y1_oxidation + initial_moles_in_stratosphere

    # After Year 2: 3% oxidation of the total remaining amount
    total_moles_at_end_of_y2 = total_moles_at_end_of_y1 * (1 - subsequent_years_oxidation_rate)

    # After Year 3: Another 3% oxidation of the total remaining amount
    final_moles_ch4_remaining = total_moles_at_end_of_y2 * (1 - subsequent_years_oxidation_rate)

    # Step 4: Calculate the total moles of air in the atmosphere
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1e3 # Convert kg to grams
    total_moles_air = total_atmosphere_mass_g / avg_molar_mass_air

    # Step 5: Calculate the final concentration increase in ppb
    # ppb = (moles of substance / total moles of mixture) * 10^9
    final_ppb_increase = (final_moles_ch4_remaining / total_moles_air) * 1e9
    
    # Step 6: Print the final equation with all its components as requested
    print("This script calculates the final increase in atmospheric methane concentration (in ppb) after 3 years.")
    print("\nFinal Equation:")
    print(f"Final Methane Moles / Total Air Moles * 1,000,000,000 = Final Concentration Increase (ppb)")
    print(f"{final_moles_ch4_remaining} / {total_moles_air} * 1000000000 = {final_ppb_increase}")
    
    # Final answer in the required format
    print("\n<<<" + f"{final_ppb_increase:.5f}" + ">>>")

calculate_methane_increase()