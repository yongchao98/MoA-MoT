import sys

def calculate_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """

    # --- 1. Define Constants ---
    total_ch4_release_tons = 250000
    total_ch4_release_kg = total_ch4_release_tons * 1000

    # Molecular weights (kg/mol)
    mw_ch4_kg = 16.0 / 1000
    mw_air_kg = 28.97 / 1000 # Standard molar mass of dry air

    # Atmospheric layer properties
    mass_troposphere_kg = 4.5e18
    mass_stratosphere_kg = 0.6e18
    temp_troposphere_k = 280
    temp_stratosphere_k = 220

    # Partitioning of the release
    tropo_share = 0.80
    strato_share = 0.20

    # Oxidation rates
    tropo_ox_y1 = 0.05
    subsequent_ox = 0.03

    # --- 2. Calculate Initial Mass Distribution ---
    tropo_initial_mass_kg = total_ch4_release_kg * tropo_share
    strato_yearly_injection_kg = (total_ch4_release_kg * strato_share) / 3

    # --- 3. Simulate Oxidation Over 3 Years ---
    # Year 1
    tropo_mass_after_y1 = tropo_initial_mass_kg * (1 - tropo_ox_y1)
    strato_mass_after_y1 = strato_yearly_injection_kg

    # Year 2
    tropo_mass_after_y2 = tropo_mass_after_y1 * (1 - subsequent_ox)
    strato_mass_after_y2 = (strato_mass_after_y1 * (1 - subsequent_ox)) + strato_yearly_injection_kg

    # Year 3
    final_tropo_mass_kg = tropo_mass_after_y2 * (1 - subsequent_ox)
    final_strato_mass_kg = (strato_mass_after_y2 * (1 - subsequent_ox)) + strato_yearly_injection_kg

    # --- 4. Calculate Final Moles in Each Layer ---
    final_ch4_moles_tropo = final_tropo_mass_kg / mw_ch4_kg
    final_ch4_moles_strato = final_strato_mass_kg / mw_ch4_kg

    # --- 5. Calculate Moles of Air in Each Layer ---
    moles_air_tropo = mass_troposphere_kg / mw_air_kg
    moles_air_strato = mass_stratosphere_kg / mw_air_kg

    # --- 6. Calculate Volume-Based Concentration (ppbv) ---
    # Numerator: Sum of (moles * temperature) for CH4
    vol_contrib_ch4 = (final_ch4_moles_tropo * temp_troposphere_k) + (final_ch4_moles_strato * temp_stratosphere_k)

    # Denominator: Sum of (moles * temperature) for Air
    vol_contrib_air = (moles_air_tropo * temp_troposphere_k) + (moles_air_strato * temp_stratosphere_k)

    # Final ppbv calculation
    ppb_increase = (vol_contrib_ch4 / vol_contrib_air) * 1e9

    # --- 7. Construct the Final Output ---
    print("Calculation of Methane Concentration Increase (ppb by volume)\n")
    print("--- Methane Mass After 3 Years ---")
    print(f"Final Methane Mass in Troposphere: {final_tropo_mass_kg:,.2f} kg")
    print(f"Final Methane Mass in Stratosphere: {final_strato_mass_kg:,.2f} kg\n")

    print("--- Moles in Each Atmospheric Layer ---")
    print(f"Final Methane Moles in Troposphere (n_CH4_tropo): {final_ch4_moles_tropo:.2e} moles")
    print(f"Final Methane Moles in Stratosphere (n_CH4_strato): {final_ch4_moles_strato:.2e} moles")
    print(f"Air Moles in Troposphere (n_air_tropo): {moles_air_tropo:.2e} moles")
    print(f"Air Moles in Stratosphere (n_air_strato): {moles_air_strato:.2e} moles\n")
    
    print("--- Final Equation for ppb Increase (Volume-Corrected) ---")
    print("ppb = ( (n_CH4_tropo * T_tropo) + (n_CH4_strato * T_strato) ) / ( (n_air_tropo * T_tropo) + (n_air_strato * T_strato) ) * 10^9\n")
    print("Plugging in the numbers:")
    print(f"ppb = ( ({final_ch4_moles_tropo:.2e} * {temp_troposphere_k}) + ({final_ch4_moles_strato:.2e} * {temp_stratosphere_k}) ) / ( ({moles_air_tropo:.2e} * {temp_troposphere_k}) + ({moles_air_strato:.2e} * {temp_stratosphere_k}) ) * 10^9")
    print(f"ppb = ( {vol_contrib_ch4:.2e} ) / ( {vol_contrib_air:.2e} ) * 10^9")
    print(f"ppb = {ppb_increase:.2f}")

    # The final answer in the specified format
    # Redirecting to stdout to be captured
    original_stdout = sys.stdout
    sys.stdout = sys.__stdout__
    print(f"\n<<<{ppb_increase:.2f}>>>")
    sys.stdout = original_stdout


calculate_methane_concentration()
