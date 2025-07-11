import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """

    # --- Given Constants ---
    total_ch4_release_tons = 250000.0  # metric tons
    mw_ch4 = 16.0  # g/mol
    mw_air = 29.0  # g/mol (average molecular weight of air)
    total_atmosphere_mass_kg = 5.1e18  # kg

    # --- Problem Conditions ---
    troposphere_fraction = 0.80
    stratosphere_fraction = 0.20
    year1_tropo_oxidation_rate = 0.05
    subsequent_oxidation_rate = 0.03
    years_for_stratospheric_mixing = 3

    # --- Step 1: Quantify Initial Methane ---
    total_ch4_release_g = total_ch4_release_tons * 1e6  # Convert metric tons to grams
    ch4_mass_to_tropo = total_ch4_release_g * troposphere_fraction
    ch4_mass_to_strato_total = total_ch4_release_g * stratosphere_fraction
    ch4_yearly_addition_to_strato = ch4_mass_to_strato_total / years_for_stratospheric_mixing

    # --- Step 2 & 3: Model Methane Mass Over 3 Years ---
    # Initialize methane mass in each layer
    ch4_in_tropo = 0.0
    ch4_in_strato = 0.0

    # Year 1
    ch4_in_tropo += ch4_mass_to_tropo
    ch4_in_strato += ch4_yearly_addition_to_strato
    # Apply oxidation for Year 1
    ch4_in_tropo *= (1 - year1_tropo_oxidation_rate)
    # No oxidation in stratosphere in year 1

    # Year 2
    ch4_in_strato += ch4_yearly_addition_to_strato
    # Apply oxidation for Year 2
    ch4_in_tropo *= (1 - subsequent_oxidation_rate)
    ch4_in_strato *= (1 - subsequent_oxidation_rate)

    # Year 3
    ch4_in_strato += ch4_yearly_addition_to_strato
    # Apply oxidation for Year 3
    ch4_in_tropo *= (1 - subsequent_oxidation_rate)
    ch4_in_strato *= (1 - subsequent_oxidation_rate)

    # --- Step 4: Calculate Total Remaining Methane ---
    total_ch4_remaining_g = ch4_in_tropo + ch4_in_strato

    # --- Step 5: Calculate Final Concentration ---
    # Convert mass to moles
    total_ch4_remaining_moles = total_ch4_remaining_g / mw_ch4
    
    # Calculate total moles of air
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1000
    total_moles_air = total_atmosphere_mass_g / mw_air

    # Calculate ppb increase
    final_ppb_increase = (total_ch4_remaining_moles / total_moles_air) * 1e9

    # --- Step 6: Output the Final Equation and Result ---
    # The format is: (moles of remaining CH4 / total moles of atmosphere) * 10^9 = final ppb
    print(f"Calculation: ({total_ch4_remaining_moles:.4g} moles CH4 / {total_moles_air:.4g} moles air) * 10^9 = {final_ppb_increase:.4f} ppb")

calculate_methane_increase()
<<<0.0804>>>