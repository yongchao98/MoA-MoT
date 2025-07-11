import sys
import io

# Redirect stdout to a string buffer to capture all prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration (in ppb)
    after 3 years of mixing and oxidation from the Nord Stream leak.
    """
    # Step 1: Define constants and initial values
    ch4_release_metric_tons = 250_000
    ch4_release_g = ch4_release_metric_tons * 1_000 * 1_000  # Convert metric tons to grams
    
    # Molecular weights (g/mol)
    mw_ch4 = 16.0
    mw_air = 29.0  # Average molecular weight of air
    
    # Atmospheric mass (kg) and conversion to grams
    mass_atm_kg = 5.1e18
    mass_atm_g = mass_atm_kg * 1_000
    
    # Oxidation rates
    oxidation_y1_tropo = 0.05
    oxidation_subsequent = 0.03
    
    # Step 2: Calculate initial moles
    total_moles_ch4_released = ch4_release_g / mw_ch4
    total_moles_air = mass_atm_g / mw_air
    
    print("Initial Calculation:")
    print(f"Total Methane Released: {total_moles_ch4_released:.4e} moles")
    print(f"Total Moles of Air in Atmosphere: {total_moles_air:.4e} moles\n")

    # Step 3: Model methane distribution and oxidation over 3 years
    
    # Initial distribution of methane moles
    moles_ch4_tropo = 0.80 * total_moles_ch4_released
    moles_ch4_strato_pool = 0.20 * total_moles_ch4_released
    moles_ch4_strato_add_per_year = moles_ch4_strato_pool / 3
    
    # Initialize stratospheric methane at zero
    moles_ch4_strato = 0.0

    print("Year-by-Year Methane Remaining (moles):")
    # Year 1
    # Oxidation only in troposphere
    moles_ch4_tropo *= (1 - oxidation_y1_tropo)
    # Mixing into stratosphere (no oxidation in strato for year 1)
    moles_ch4_strato += moles_ch4_strato_add_per_year
    print(f"End of Year 1: Troposphere={moles_ch4_tropo:.4e}, Stratosphere={moles_ch4_strato:.4e}")

    # Year 2
    # Oxidation in troposphere
    moles_ch4_tropo *= (1 - oxidation_subsequent)
    # Mixing into stratosphere, then oxidation
    moles_ch4_strato += moles_ch4_strato_add_per_year
    moles_ch4_strato *= (1 - oxidation_subsequent)
    print(f"End of Year 2: Troposphere={moles_ch4_tropo:.4e}, Stratosphere={moles_ch4_strato:.4e}")
    
    # Year 3
    # Oxidation in troposphere
    moles_ch4_tropo *= (1 - oxidation_subsequent)
    # Mixing into stratosphere, then oxidation
    moles_ch4_strato += moles_ch4_strato_add_per_year
    moles_ch4_strato *= (1 - oxidation_subsequent)
    print(f"End of Year 3: Troposphere={moles_ch4_tropo:.4e}, Stratosphere={moles_ch4_strato:.4e}\n")

    # Step 4: Calculate final concentration
    final_moles_ch4 = moles_ch4_tropo + moles_ch4_strato
    ppb_increase = (final_moles_ch4 / total_mles_air) * 1e9

    # Step 5: Format the output
    print("Final Calculation:")
    print(f"Total remaining methane after 3 years: {final_moles_ch4:.4e} moles")
    print(f"Total moles of air: {total_moles_air:.4e} moles")
    print("Equation: (Final Methane Moles / Total Air Moles) * 1,000,000,000")
    print(f"Result: ({final_moles_ch4:.4e} / {total_moles_air:.4e}) * 1e9 = {ppb_increase:.4f} ppb")
    
    return ppb_increase

# Execute the calculation
final_ppb = calculate_methane_increase()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# Final answer in the required format
print(f"<<<{final_ppb:.4f}>>>")