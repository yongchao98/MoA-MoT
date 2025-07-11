import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """

    # --- Given Constants ---
    # Mass of CH4 released
    mass_ch4_released_tons = 250000  # metric tons
    # Total mass of the atmosphere
    mass_atmosphere_kg = 5.1e18  # kg
    # Molar mass of CH4
    molar_mass_ch4_g_mol = 16.0  # g/mol
    # Average molar mass of air (approximated)
    molar_mass_air_g_mol = 29.0  # g/mol
    # Oxidation reduction rates
    reduction_year_1 = 0.05  # 5%
    reduction_subsequent_years = 0.03  # 3%

    # --- Step 1: Calculate the initial increase in ppb ---

    # Convert masses to consistent units (grams)
    mass_ch4_released_g = mass_ch4_released_tons * 1e6  # 1 metric ton = 1e6 grams
    mass_atmosphere_g = mass_atmosphere_kg * 1e3     # 1 kg = 1e3 grams

    # Calculate moles of CH4 released
    moles_ch4 = mass_ch4_released_g / molar_mass_ch4_g_mol

    # Calculate total moles of air in the atmosphere
    moles_atmosphere = mass_atmosphere_g / molar_mass_air_g_mol

    # Calculate the initial increase in mole fraction and convert to ppb
    initial_increase_ppb = (moles_ch4 / moles_atmosphere) * 1e9

    print(f"Calculation for Initial Methane Increase:")
    print(f"Moles of CH4 Released = {mass_ch4_released_g:.2e} g / {molar_mass_ch4_g_mol} g/mol = {moles_ch4:.4e} mol")
    print(f"Total Moles of Atmosphere = {mass_atmosphere_g:.2e} g / {molar_mass_air_g_mol} g/mol = {moles_atmosphere:.4e} mol")
    print(f"Initial Increase (ppb) = ({moles_ch4:.4e} mol / {moles_atmosphere:.4e} mol) * 1e9 = {initial_increase_ppb:.4f} ppb\n")

    # --- Step 2: Apply annual oxidation reductions ---

    print("Calculating concentration after annual oxidation:")
    # Concentration after 1 year
    ppb_after_1_year = initial_increase_ppb * (1 - reduction_year_1)
    print(f"After Year 1 (5% reduction): {initial_increase_ppb:.4f} * (1 - {reduction_year_1}) = {ppb_after_1_year:.4f} ppb")

    # Concentration after 2 years
    ppb_after_2_years = ppb_after_1_year * (1 - reduction_subsequent_years)
    print(f"After Year 2 (3% reduction): {ppb_after_1_year:.4f} * (1 - {reduction_subsequent_years}) = {ppb_after_2_years:.4f} ppb")

    # Concentration after 3 years
    ppb_after_3_years = ppb_after_2_years * (1 - reduction_subsequent_years)
    print(f"After Year 3 (3% reduction): {ppb_after_2_years:.4f} * (1 - {reduction_subsequent_years}) = {ppb_after_3_years:.4f} ppb\n")

    # --- Final Answer ---
    print("Final Answer:")
    print(f"The final increase in atmospheric methane concentration after 3 years is {ppb_after_3_years:.4f} ppb.")
    
    # The final answer is returned in the specified format.
    # Using round to avoid floating point inaccuracies in the final output string.
    final_answer = round(ppb_after_3_years, 4)
    return f"<<<{final_answer}>>>"

# Execute the function and print the final formatted answer
final_result_string = calculate_methane_increase()
print(final_result_string)