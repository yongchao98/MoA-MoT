import sys

def solve_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration after 3 years
    based on the Nord Stream leak data and subsequent atmospheric processes.
    """

    # --- Step 0: Define initial constants from the problem ---
    total_ch4_release_tonnes = 250000
    total_ch4_release_kg = total_ch4_release_tonnes * 1000

    # Mixing parameters
    tropo_mix_fraction = 0.80
    strato_mix_fraction = 0.20
    strato_mixing_period_years = 3

    # Oxidation parameters
    oxidation_y1_tropo = 0.05
    oxidation_y2_y3_total = 0.03

    # Atmospheric properties
    total_atm_mass_kg = 5.1e18
    mw_air_g_mol = 29.0  # Average molecular weight of air
    mw_ch4_g_mol = 16.0  # Molecular weight of methane

    # --- Step 1: Calculate annual mass distribution ---
    ch4_to_tropo_initial = total_ch4_release_kg * tropo_mix_fraction
    ch4_to_strato_per_year = (total_ch4_release_kg * strato_mix_fraction) / strato_mixing_period_years

    # --- Step 2: Model the mass change year-by-year ---

    # Year 1
    mass_tropo_after_y1_oxidation = ch4_to_tropo_initial * (1 - oxidation_y1_tropo)
    mass_end_y1 = mass_tropo_after_y1_oxidation + ch4_to_strato_per_year

    # Year 2
    mass_before_y2_oxidation = mass_end_y1 + ch4_to_strato_per_year
    mass_end_y2 = mass_before_y2_oxidation * (1 - oxidation_y2_y3_total)

    # Year 3
    mass_before_y3_oxidation = mass_end_y2 + ch4_to_strato_per_year
    final_ch4_mass_kg = mass_before_y3_oxidation * (1 - oxidation_y2_y3_total)

    # --- Step 3: Calculate the final concentration increase in ppb ---
    mass_fraction = final_ch4_mass_kg / total_atm_mass_kg
    mole_fraction_ratio = mw_air_g_mol / mw_ch4_g_mol
    final_ppb_increase = mass_fraction * mole_fraction_ratio * 1e9

    # --- Step 4: Print the step-by-step calculation ---
    print("Step 1: Calculate methane mass at the end of Year 1")
    print(f"Mass after 1 year = (Initial Tropospheric Mass * (1 - Tropospheric Oxidation)) + Stratospheric Injection")
    print(f"Mass after 1 year = ({ch4_to_tropo_initial:.2e} kg * (1 - {oxidation_y1_tropo})) + {ch4_to_strato_per_year:.2e} kg = {mass_end_y1:.2e} kg\n")

    print("Step 2: Calculate methane mass at the end of Year 2")
    print(f"Mass after 2 years = (Mass from Year 1 + Stratospheric Injection) * (1 - Total Oxidation)")
    print(f"Mass after 2 years = ({mass_end_y1:.2e} kg + {ch4_to_strato_per_year:.2e} kg) * (1 - {oxidation_y2_y3_total}) = {mass_end_y2:.2e} kg\n")

    print("Step 3: Calculate methane mass at the end of Year 3")
    print(f"Mass after 3 years = (Mass from Year 2 + Stratospheric Injection) * (1 - Total Oxidation)")
    print(f"Mass after 3 years = ({mass_end_y2:.2e} kg + {ch4_to_strato_per_year:.2e} kg) * (1 - {oxidation_y2_y3_total}) = {final_ch4_mass_kg:.2e} kg\n")

    print("Step 4: Calculate final concentration increase in ppb")
    print(f"ppb Increase = (Final CH4 Mass / Total Atmosphere Mass) * (MW_Air / MW_CH4) * 1,000,000,000")
    print(f"ppb Increase = ({final_ch4_mass_kg:.2e} kg / {total_atm_mass_kg:.2e} kg) * ({mw_air_g_mol} / {mw_ch4_g_mol}) * 1e9")
    print(f"ppb Increase = {final_ppb_increase:.6f}\n")
    
    # Use this format to avoid scientific notation for the final answer to be extracted.
    sys.stdout.write(f'<<<{final_ppb_increase:.6f}>>>')

solve_methane_concentration()