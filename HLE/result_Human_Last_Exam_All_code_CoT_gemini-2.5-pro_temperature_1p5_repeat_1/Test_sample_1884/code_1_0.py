def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years
    following the Nord Stream leak, accounting for mixing and oxidation.
    """

    # --- Step 1: Define initial conditions and constants ---
    total_release_tons = 250000
    total_release_kg = total_release_tons * 1000  # Convert metric tons to kg

    # Masses and Molar Masses
    total_atmosphere_mass_kg = 5.1e18
    molar_mass_ch4 = 16.0  # g/mol
    molar_mass_air = 29.0   # g/mol, a standard approximation for dry air

    # Mixing and Oxidation Parameters
    tropospheric_fraction = 0.80
    stratospheric_fraction = 0.20
    oxidation_year_1 = 0.05
    oxidation_subsequent_years = 0.03

    # --- Step 2: Calculate the mass distribution ---
    mass_tropo_initial = total_release_kg * tropospheric_fraction
    mass_strato_pool = total_release_kg * stratospheric_fraction
    
    # --- Step 3: Track the mass remaining from the tropospheric portion ---
    # Apply decay for Year 1, Year 2, and Year 3
    mass_tropo_final = mass_tropo_initial * (1 - oxidation_year_1) * (1 - oxidation_subsequent_years) * (1 - oxidation_subsequent_years)

    # --- Step 4: Track the mass remaining from the stratospheric portion year by year ---
    strato_add_per_year = mass_strato_pool / 3.0
    
    # End of Year 1: Addition without oxidation
    mass_strato_y1 = strato_add_per_year
    
    # End of Year 2: Add more mass, then apply oxidation to the total
    mass_strato_y2 = (mass_strato_y1 + strato_add_per_year) * (1 - oxidation_subsequent_years)
    
    # End of Year 3: Add the final portion, then apply oxidation to the new total
    mass_strato_final = (mass_strato_y2 + strato_add_per_year) * (1 - oxidation_subsequent_years)
    
    # --- Step 5: Sum the final masses ---
    total_remaining_mass_kg = mass_tropo_final + mass_strato_final
    
    # --- Step 6: Convert the final mass to ppb increase ---
    # The formula is: ppb = (mass_fraction) * (molar_mass_ratio) * 1e9
    ppb_increase = (total_remaining_mass_kg / total_atmosphere_mass_kg) * (molar_mass_air / molar_mass_ch4) * 1e9

    # --- Step 7: Print the results ---
    print("Calculation of Methane Concentration Increase After 3 Years:\n")
    print(f"Initial methane release: {total_release_kg:,.0f} kg")
    print(f"Remaining mass after 3 years: {total_remaining_mass_kg:,.2f} kg\n")
    
    print("Final concentration is calculated using the formula:")
    print("ppb = (Total Remaining CH4 Mass / Total Atmosphere Mass) * (Molar Mass Air / Molar Mass CH4) * 1,000,000,000\n")
    
    print("Values for the final equation:")
    print(f"Total Remaining CH4 Mass = {total_remaining_mass_kg:.2f} kg")
    print(f"Total Atmosphere Mass = {total_atmosphere_mass_kg:.1e} kg")
    print(f"Molar Mass Air = {molar_mass_air} g/mol")
    print(f"Molar Mass CH4 = {molar_mass_ch4} g/mol")
    
    print("\nResult:")
    print(f"({total_remaining_mass_kg:.2f} / {total_atmosphere_mass_kg:.1e}) * ({molar_mass_air} / {molar_mass_ch4}) * 1,000,000,000 = {ppb_increase:.5f} ppb")
    
    # The final answer in the required format
    print(f"\nFinal calculated increase in atmospheric methane concentration: {ppb_increase:.5f} ppb")
    
    return ppb_increase

# Execute the calculation
final_ppb_increase = calculate_methane_increase()
# The final answer is wrapped in <<<>>>
print(f"\n<<<>>>\n{final_ppb_increase:.5f}")