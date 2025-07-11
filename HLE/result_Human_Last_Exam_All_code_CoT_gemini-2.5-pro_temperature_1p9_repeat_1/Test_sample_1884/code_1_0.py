import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """
    # 1. Define Constants
    # Given values from the problem description
    total_release_tons = 250000  # metric tons
    molar_mass_ch4 = 16.0  # g/mol
    molar_mass_air = 29.0  # Approximate average molar mass of air, g/mol
    mass_troposphere_kg = 4.5e18  # kg
    mass_stratosphere_kg = 0.6e18 # kg
    total_atmosphere_mass_kg = 5.1e18 # kg

    # Mixing and oxidation parameters
    fraction_to_troposphere = 0.80
    oxidation_year1_troposphere = 0.05
    oxidation_subsequent_years = 0.03

    # Unit conversions
    total_release_g = total_release_tons * 1e6  # tons to grams
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1e3 # kg to grams

    # 2. Initial Mass Distribution
    initial_mass_troposphere_g = total_release_g * fraction_to_troposphere
    initial_mass_stratosphere_g = total_release_g * (1 - fraction_to_troposphere)

    # 3. Calculate Methane Decay over 3 years
    # Tropospheric mass decay: 1 year at 5% loss, 2 years at 3% loss
    final_mass_troposphere_g = initial_mass_troposphere_g * (1 - oxidation_year1_troposphere) * ((1 - oxidation_subsequent_years) ** 2)

    # Stratospheric mass decay: 2 years at 3% loss (for years 2 and 3)
    final_mass_stratosphere_g = initial_mass_stratosphere_g * ((1 - oxidation_subsequent_years) ** 2)

    # 4. Calculate Total Remaining Methane
    total_final_mass_g = final_mass_troposphere_g + final_mass_stratosphere_g

    # 5. Calculate Final Concentration Increase in ppb
    # ppb is a mole fraction, so we convert from mass fraction using molar masses.
    # ppb_increase = (mass_CH4 / mass_Air) * (M_Air / M_CH4) * 10^9
    final_ppb_increase = (total_final_mass_g / total_atmosphere_mass_g) * (molar_mass_air / molar_mass_ch4) * 1e9

    # 6. Display the Final Equation and Result
    # This fulfills the requirement to output each number in the final equation.
    print("Calculation for the final increase in atmospheric CH4 concentration (ppb):")
    print(f"Final CH4 Mass = {total_final_mass_g:.4e} g")
    print(f"Total Atmosphere Mass = {total_atmosphere_mass_g:.2e} g")
    print(f"Molar Mass of Air = {molar_mass_air} g/mol")
    print(f"Molar Mass of CH4 = {molar_mass_ch4} g/mol")
    print("\nFinal Equation:")
    print(f"Increase (ppb) = ({total_final_mass_g:.4e} g / {total_atmosphere_mass_g:.2e} g) * ({molar_mass_air} / {molar_mass_ch4}) * 1,000,000,000")
    print(f"\nCalculated Increase: {final_ppb_increase:.4f} ppb")
    
    # Final answer in the specified format
    print(f"<<<{final_ppb_increase:.4f}>>>")

solve_methane_concentration()

# Restore original stdout and print captured output
sys.stdout = original_stdout
output_str = captured_output.getvalue()
print(output_str)
