import sys

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years
    based on the Nord Stream pipeline leak scenario.
    """
    # --- Step 1: Define Constants and Initial Values ---

    total_methane_release_tons = 250000
    total_methane_release_kg = total_methane_release_tons * 1000

    # Split the release into two pulses
    tropospheric_pulse_kg = total_methane_release_kg * 0.80
    gradual_pulse_kg = total_methane_release_kg * 0.20
    gradual_addition_per_year_kg = gradual_pulse_kg / 3

    # Atmospheric and chemical constants
    molar_mass_ch4_kg_per_mol = 0.016  # 16 g/mol
    molar_mass_air_kg_per_mol = 0.029  # ~29 g/mol, standard approximation
    total_mass_atmosphere_kg = 5.1e18

    # Oxidation rates
    oxidation_rate_y1 = 0.05  # 5%
    oxidation_rate_y2_y3 = 0.03  # 3%

    print("--- Calculation Plan ---")
    print(f"Total Methane Released: {total_methane_release_kg:,.0f} kg")
    print(f"Initial Tropospheric Pulse (80%): {tropospheric_pulse_kg:,.0f} kg")
    print(f"Gradual Addition per year (20%/3): {gradual_addition_per_year_kg:,.2f} kg\n")

    # --- Step 2: Track Methane Mass Year-by-Year ---

    # Year 1
    # The 80% pulse is oxidized by 5% in the troposphere. The first gradual addition occurs.
    methane_mass_after_y1_ox = tropospheric_pulse_kg * (1 - oxidation_rate_y1)
    total_mass_end_y1 = methane_mass_after_y1_ox + gradual_addition_per_year_kg

    print("--- Year 1 Calculation ---")
    print(f"Equation: (Tropospheric Pulse * (1 - Oxidation Y1)) + Gradual Addition")
    print(f"Result: ({tropospheric_pulse_kg:,.0f} kg * (1 - {oxidation_rate_y1})) + {gradual_addition_per_year_kg:,.2f} kg = {total_mass_end_y1:,.2f} kg remaining\n")

    # Year 2
    # The total from Y1 and the second gradual addition are oxidized by 3%.
    mass_before_y2_ox = total_mass_end_y1 + gradual_addition_per_year_kg
    total_mass_end_y2 = mass_before_y2_ox * (1 - oxidation_rate_y2_y3)

    print("--- Year 2 Calculation ---")
    print(f"Equation: (Mass from Y1 + Gradual Addition) * (1 - Oxidation Y2)")
    print(f"Result: ({total_mass_end_y1:,.2f} kg + {gradual_addition_per_year_kg:,.2f} kg) * (1 - {oxidation_rate_y2_y3}) = {total_mass_end_y2:,.2f} kg remaining\n")

    # Year 3
    # The total from Y2 and the third gradual addition are oxidized by 3%.
    mass_before_y3_ox = total_mass_end_y2 + gradual_addition_per_year_kg
    final_methane_mass_kg = mass_before_y3_ox * (1 - oxidation_rate_y2_y3)

    print("--- Year 3 Calculation ---")
    print(f"Equation: (Mass from Y2 + Gradual Addition) * (1 - Oxidation Y3)")
    print(f"Result: ({total_mass_end_y2:,.2f} kg + {gradual_addition_per_year_kg:,.2f} kg) * (1 - {oxidation_rate_y2_y3}) = {final_methane_mass_kg:,.2f} kg remaining\n")


    # --- Step 3: Calculate Final Concentration in ppb ---
    
    # Convert final mass of CH4 to moles
    final_moles_ch4 = final_methane_mass_kg / molar_mass_ch4_kg_per_mol

    # Calculate total moles of air
    total_moles_air = total_mass_atmosphere_kg / molar_mass_air_kg_per_mol

    # Calculate final ppb increase
    ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9

    print("--- Final Concentration Calculation ---")
    print("Equation: (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
    print(f"Final Moles CH4 = {final_methane_mass_kg:,.2f} kg / {molar_mass_ch4_kg_per_mol} kg/mol = {final_moles_ch4:,.2f} moles")
    print(f"Total Moles Air = {total_mass_atmosphere_kg:.1e} kg / {molar_mass_air_kg_per_mol} kg/mol = {total_moles_air:.3e} moles")
    print(f"Final ppb Increase = ({final_moles_ch4:.2f} / {total_moles_air:.3e}) * 1e9 = {ppb_increase:.2f} ppb\n")
    
    # Final answer in the required format
    sys.stdout.write(f"<<<{ppb_increase:.2f}>>>")

calculate_methane_increase()