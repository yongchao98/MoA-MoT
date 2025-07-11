import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years
    following the Nord Stream pipeline leaks.
    """

    # --- Step 1: Constants and Initial Values ---
    mass_ch4_released_kg = 250000 * 1000  # 250,000 metric tons to kg
    molar_mass_ch4_kg_mol = 16.0 / 1000   # g/mol to kg/mol
    total_mass_atmosphere_kg = 5.1e18
    # Assumption: Average molar mass of air is ~29 g/mol
    molar_mass_air_kg_mol = 29.0 / 1000

    # --- Step 2: Calculate total moles of CH4 released ---
    moles_released = mass_ch4_released_kg / molar_mass_ch4_kg_mol

    # --- Step 3: Model Tropospheric Methane ---
    moles_tropo_initial = 0.80 * moles_released
    # Year 1: 5% reduction
    moles_tropo_after_y1 = moles_tropo_initial * (1 - 0.05)
    # Year 2: 3% reduction
    moles_tropo_after_y2 = moles_tropo_after_y1 * (1 - 0.03)
    # Year 3: 3% reduction
    moles_tropo_final = moles_tropo_after_y2 * (1 - 0.03)

    # --- Step 4: Model Stratospheric Methane ---
    moles_strato_pool = 0.20 * moles_released
    strato_annual_influx = moles_strato_pool / 3.0
    
    # End of Year 1: Influx occurs, but "subsequent year" oxidation hasn't started
    moles_strato_after_y1 = strato_annual_influx
    
    # End of Year 2: Amount from Y1 + new influx, then 3% oxidation
    moles_strato_after_y2 = (moles_strato_after_y1 + strato_annual_influx) * (1 - 0.03)
    
    # End of Year 3: Amount from Y2 + new influx, then 3% oxidation
    moles_strato_final = (moles_strato_after_y2 + strato_annual_influx) * (1 - 0.03)

    # --- Step 5: Sum Remaining Methane ---
    total_moles_remaining = moles_tropo_final + moles_strato_final

    # --- Step 6: Calculate Total Air Moles ---
    total_moles_air = total_mass_atmosphere_kg / molar_mass_air_kg_mol

    # --- Step 7: Calculate Final Concentration Increase ---
    increase_ppb = (total_moles_remaining / total_moles_air) * 1e9

    # --- Final Output ---
    print("After 3 years of mixing and oxidation:")
    print(f"Remaining methane from leak in troposphere: {moles_tropo_final:.4e} moles")
    print(f"Remaining methane from leak in stratosphere: {moles_strato_final:.4e} moles")
    print("-" * 40)
    print(f"Total remaining moles of leaked CH4 = {total_moles_remaining:.4f} moles")
    print(f"Total moles of air in atmosphere = {total_moles_air:.4f} moles")
    print("\nFinal Calculation for the concentration increase:")
    print(f"({total_moles_remaining:.4f} moles_CH4 / {total_moles_air:.4f} moles_air) * 1,000,000,000 = {increase_ppb:.3f} ppb")


calculate_methane_increase()
<<<80.425>>>