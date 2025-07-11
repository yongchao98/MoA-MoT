import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """

    # --- Step 1: Define Initial Values ---
    total_ch4_release_tons = 250000
    total_ch4_release_kg = float(total_ch4_release_tons * 1000)

    molar_mass_ch4 = 16.0  # g/mol
    molar_mass_air = 28.97 # g/mol (average)

    mass_atmosphere_total = 5.1e18 # kg

    # Oxidation rates
    oxidation_rate_y1_tropo = 0.05
    oxidation_rate_subsequent = 0.03

    # --- Step 2: Calculate Methane Distribution ---
    ch4_in_tropo_initial = total_ch4_release_kg * 0.80
    ch4_to_strato_total = total_ch4_release_kg * 0.20
    ch4_to_strato_per_year = ch4_to_strato_total / 3

    # --- Step 3: Simulate Year-by-Year Changes ---
    ch4_in_tropo_current = 0.0
    ch4_in_strato_current = 0.0

    # --- Year 1 ---
    # Methane added to troposphere, then oxidized
    ch4_in_tropo_current = ch4_in_tropo_initial * (1 - oxidation_rate_y1_tropo)
    # Methane added to stratosphere (no oxidation in year 1 for this part)
    ch4_in_strato_current += ch4_to_strato_per_year

    # --- Year 2 ---
    # Existing tropospheric methane is oxidized
    ch4_in_tropo_current *= (1 - oxidation_rate_subsequent)
    # Methane added to stratosphere, then total stratospheric amount is oxidized
    ch4_in_strato_current += ch4_to_strato_per_year
    ch4_in_strato_current *= (1 - oxidation_rate_subsequent)

    # --- Year 3 ---
    # Existing tropospheric methane is oxidized
    ch4_in_tropo_current *= (1 - oxidation_rate_subsequent)
    # Methane added to stratosphere, then total stratospheric amount is oxidized
    ch4_in_strato_current += ch4_to_strato_per_year
    ch4_in_strato_current *= (1 - oxidation_rate_subsequent)

    # --- Step 4: Calculate Final Total Methane Mass ---
    total_ch4_final_kg = ch4_in_tropo_current + ch4_in_strato_current

    # --- Step 5 & 6: Convert to ppb and Output ---
    final_mass_fraction = total_ch4_final_kg / mass_atmosphere_total
    ppb_increase = final_mass_fraction * (molar_mass_air / molar_mass_ch4) * 1e9
    
    print("--- Final Calculation ---")
    print(f"Remaining methane mass after 3 years: {total_ch4_final_kg:.2f} kg")
    print(f"Total mass of the atmosphere: {mass_atmosphere_total:.1e} kg\n")
    
    print("The final concentration increase is calculated as:")
    print("ppb = (Total Remaining CH₄ Mass / Total Atmosphere Mass) * (Molar Mass of Air / Molar Mass of CH₄) * 1,000,000,000\n")
    
    print("Plugging in the numbers:")
    # Using 'g' format for cleaner representation of scientific notation
    print(f"ppb = ({total_ch4_final_kg:.4g} kg / {mass_atmosphere_total:.4g} kg) * ({molar_mass_air} g/mol / {molar_mass_ch4} g/mol) * 1e9")
    
    print(f"\nFinal Increase in Atmospheric Methane Concentration: {ppb_increase:.3f} ppb")
    
    # The final answer for the system
    print(f"\n<<<{ppb_increase:.3f}>>>")

calculate_methane_increase()