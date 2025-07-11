import math

def calculate_methane_concentration_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years
    of mixing and oxidation from the Nord Stream pipeline leak.
    """

    # --- Step 1: Define constants and initial values ---

    # Given values
    methane_release_metric_tons = 250000
    mw_ch4_g_per_mol = 16.0
    total_atmosphere_mass_kg = 5.1e18
    
    # Assumed standard value
    avg_mw_air_g_per_mol = 29.0

    # Conversion factors
    g_per_kg = 1000
    kg_per_metric_ton = 1000

    # Oxidation rates
    oxidation_rate_y1_troposphere = 0.05  # 5%
    oxidation_rate_subsequent_years = 0.03  # 3%

    # --- Step 2: Calculate initial moles ---

    # Convert released methane mass to moles
    methane_release_g = methane_release_metric_tons * kg_per_metric_ton * g_per_kg
    initial_moles_ch4 = methane_release_g / mw_ch4_g_per_mol

    # Convert total atmosphere mass to moles
    total_atmosphere_mass_g = total_atmosphere_mass_kg * g_per_kg
    total_moles_atmosphere = total_atmosphere_mass_g / avg_mw_air_g_per_mol
    
    # Initial distribution of the leaked methane
    moles_for_troposphere = initial_moles_ch4 * 0.80
    moles_for_stratosphere_pool = initial_moles_ch4 * 0.20
    
    # Stratospheric mixing occurs over 3 years
    stratospheric_mixing_per_year = moles_for_stratosphere_pool / 3

    # --- Step 3: Simulate methane dynamics over 3 years ---

    # Initialize the state variables for methane moles in each reservoir
    moles_in_troposphere = 0
    moles_in_stratosphere = 0
    moles_in_pool = moles_for_stratosphere_pool

    # Year 1 Simulation
    # Mixing
    moles_in_troposphere += moles_for_troposphere
    moles_mixed_to_stratosphere_y1 = stratospheric_mixing_per_year
    moles_in_stratosphere += moles_mixed_to_stratosphere_y1
    moles_in_pool -= moles_mixed_to_stratosphere_y1
    # Oxidation (only in troposphere for Y1)
    moles_in_troposphere *= (1 - oxidation_rate_y1_troposphere)
    
    # Year 2 Simulation
    # Mixing
    moles_mixed_to_stratosphere_y2 = stratospheric_mixing_per_year
    moles_in_stratosphere += moles_mixed_to_stratosphere_y2
    moles_in_pool -= moles_mixed_to_stratosphere_y2
    # Oxidation (in both troposphere and stratosphere)
    moles_in_troposphere *= (1 - oxidation_rate_subsequent_years)
    moles_in_stratosphere *= (1 - oxidation_rate_subsequent_years)

    # Year 3 Simulation
    # Mixing
    moles_mixed_to_stratosphere_y3 = stratospheric_mixing_per_year
    moles_in_stratosphere += moles_mixed_to_stratosphere_y3
    mles_in_pool = 0 # Pool is now empty
    # Oxidation (in both troposphere and stratosphere)
    moles_in_troposphere *= (1 - oxidation_rate_subsequent_years)
    moles_in_stratosphere *= (1 - oxidation_rate_subsequent_years)
    
    # --- Step 4: Calculate final concentration ---
    
    final_moles_ch4 = moles_in_troposphere + moles_in_stratosphere
    
    # Final concentration increase in ppb (parts per billion)
    final_ppb_increase = (final_moles_ch4 / total_moles_atmosphere) * 1e9

    # --- Step 5: Format the output ---

    print("Calculation Breakdown:")
    print(f"Initial moles of CH₄ released: {initial_moles_ch4:.4g}")
    print(f"Total moles of the atmosphere: {total_moles_atmosphere:.4g}")
    print(f"Final moles of CH₄ remaining after 3 years: {final_moles_ch4:.4g}")
    print("\nFinal Equation for Concentration Increase (ppb):")
    print(f"(Final Moles CH₄ / Total Moles Atmosphere) * 1,000,000,000")
    print(f"= ({final_moles_ch4:.10f} / {total_moles_atmosphere:.10f}) * 1,000,000,000")
    print(f"= {final_ppb_increase:.6f} ppb")
    
    # Return final answer in the specified format
    print(f"\n<<<{final_ppb_increase:.4f}>>>")

calculate_methane_concentration_increase()