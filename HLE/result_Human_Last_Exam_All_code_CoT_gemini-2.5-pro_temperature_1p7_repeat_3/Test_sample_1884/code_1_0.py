import sys

def solve_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years
    following the Nord Stream leak, considering distribution and oxidation.
    """

    # --- Step 1: Define Constants and Initial Values ---
    total_ch4_release_tons = 250000.0
    ch4_molecular_weight = 16.0  # g/mol
    total_atmosphere_mass_kg = 5.1e18
    # Using the standard value for the average molecular weight of dry air, as it's necessary for the calculation.
    air_molecular_weight = 28.97  # g/mol

    # --- Step 2: Calculate Methane Mass Remaining After 3 Years ---

    # Methane distribution plan
    tropo_share_tons = total_ch4_release_tons * 0.80
    strato_share_total_tons = total_ch4_release_tons * 0.20
    strato_per_year_tons = strato_share_total_tons / 3.0

    print("--- Methane Mass Calculation (in metric tons) ---")

    # Year 1: 80% mixes into troposphere (with 5% oxidation), 1/3 of the 20% mixes into stratosphere
    mass_after_y1 = (tropo_share_tons * (1 - 0.05)) + strato_per_year_tons
    print(f"End of Year 1 Mass: ({tropo_share_tons:.2f} * 0.95) + {strato_per_year_tons:.2f} = {mass_after_y1:.2f}")

    # Year 2: Carry over mass, add another 1/3 of stratospheric share, then 3% total oxidation
    mass_before_y2_ox = mass_after_y1 + strato_per_year_tons
    mass_after_y2 = mass_before_y2_ox * (1 - 0.03)
    print(f"End of Year 2 Mass: ({mass_after_y1:.2f} + {strato_per_year_tons:.2f}) * 0.97 = {mass_after_y2:.2f}")

    # Year 3: Carry over mass, add final 1/3 of stratospheric share, then 3% total oxidation
    mass_before_y3_ox = mass_after_y2 + strato_per_year_tons
    final_ch4_mass_tons = mass_before_y3_ox * (1 - 0.03)
    print(f"End of Year 3 Mass: ({mass_after_y2:.2f} + {strato_per_year_tons:.2f}) * 0.97 = {final_ch4_mass_tons:.2f}")
    print("-" * 50)

    # --- Step 3: Convert Masses to Moles ---
    # Convert final methane mass from tons -> grams -> moles
    final_ch4_mass_g = final_ch4_mass_tons * 1e6
    final_ch4_moles = final_ch4_mass_g / ch4_molecular_weight

    # Convert total atmosphere mass from kg -> grams -> moles
    atmosphere_mass_g = total_atmosphere_mass_kg * 1e3
    total_moles_atmosphere = atmosphere_mass_g / air_molecular_weight

    print("--- Moles Calculation ---")
    print(f"Final Moles of CH4: ({final_ch4_mass_tons:.2f} tons * 1,000,000 g/ton) / {ch4_molecular_weight} g/mol = {final_ch4_moles:.2e} moles")
    print(f"Total Moles of Atmosphere: ({total_atmosphere_mass_kg:.2e} kg * 1,000 g/kg) / {air_molecular_weight} g/mol = {total_moles_atmosphere:.2e} moles")
    print("-" * 50)

    # --- Step 4: Calculate Final Concentration Increase in ppb ---
    increase_in_ppb = (final_ch4_moles / total_moles_atmosphere) * 1e9

    print("--- Final ppb Increase Calculation ---")
    print(f"Equation: ({final_ch4_moles:.2e} moles CH4 / {total_moles_atmosphere:.2e} total moles) * 1,000,000,000")
    print(f"Result: {increase_in_ppb:.4f} ppb")
    
    # --- Step 5: Final Answer ---
    # The 'file=sys.stderr' argument is used to separate the final answer from the explanatory text.
    print(f"<<<{increase_in_ppb:.4f}>>>", file=sys.stderr)

solve_methane_concentration()