def calculate_methane_concentration_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years,
    based on the Nord Stream leak event and subsequent atmospheric processes.
    """

    # --- Given Constants ---
    initial_ch4_release_tons = 250000  # metric tons
    total_atmosphere_mass_kg = 5.1e18 # kg
    mw_ch4 = 16.0                     # g/mol
    mw_air = 29.0                     # g/mol (average molecular weight of air)

    # --- Step 1: Calculate the initial mass of methane released in grams ---
    initial_ch4_mass_g = initial_ch4_release_tons * 1_000_000
    print(f"Step 1: Initial Methane Mass Calculation")
    print(f"The total released methane is {initial_ch4_release_tons:,} metric tons, which is equal to {initial_ch4_mass_g:e} g.\n")

    # --- Step 2: Calculate methane mass after Year 1 mixing and oxidation ---
    tropospheric_ch4_mass_g = initial_ch4_mass_g * 0.80
    stratospheric_ch4_unmixed_g = initial_ch4_mass_g * 0.20
    
    # 5% oxidation in the troposphere
    tropospheric_ch4_after_oxidation_g = tropospheric_ch4_mass_g * (1 - 0.05)
    
    total_ch4_end_y1_g = tropospheric_ch4_after_oxidation_g + stratospheric_ch4_unmixed_g
    print(f"Step 2: Methane Mass after Year 1")
    print(f"Tropospheric portion (80%): {tropospheric_ch4_mass_g:e} g undergoes a 5% reduction to {tropospheric_ch4_after_oxidation_g:e} g.")
    print(f"The unmixed stratospheric portion (20%) remains: {stratospheric_ch4_unmixed_g:e} g.")
    print(f"Total methane at end of Year 1 = {tropospheric_ch4_after_oxidation_g:.2e} g + {stratospheric_ch4_unmixed_g:.2e} g = {total_ch4_end_y1_g:e} g.\n")
    
    # --- Step 3: Calculate methane mass after Year 2 oxidation ---
    # 3% oxidation on the total remaining methane
    total_ch4_end_y2_g = total_ch4_end_y1_g * (1 - 0.03)
    print(f"Step 3: Methane Mass after Year 2")
    print(f"The mass at the start of Year 2 ({total_ch4_end_y1_g:e} g) is reduced by 3%.")
    print(f"Total methane at end of Year 2 = {total_ch4_end_y1_g:.2e} g * 0.97 = {total_ch4_end_y2_g:e} g.\n")
    
    # --- Step 4: Calculate methane mass after Year 3 oxidation ---
    # Another 3% oxidation on the total remaining methane
    final_ch4_mass_g = total_ch4_end_y2_g * (1 - 0.03)
    print(f"Step 4: Methane Mass after Year 3")
    print(f"The mass at the start of Year 3 ({total_ch4_end_y2_g:e} g) is reduced by another 3%.")
    print(f"Final remaining methane mass = {total_ch4_end_y2_g:.2e} g * 0.97 = {final_ch4_mass_g:e} g.\n")
    
    # --- Step 5: Calculate the final concentration increase in ppb ---
    # Convert final methane mass to moles
    final_ch4_moles = final_ch4_mass_g / mw_ch4
    
    # Calculate total moles in the atmosphere
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1000
    total_atmosphere_moles = total_atmosphere_mass_g / mw_air
    
    # Calculate concentration increase in ppb (parts per billion)
    ppb_increase = (final_ch4_moles / total_atmosphere_moles) * 1e9
    
    print("Step 5: Final Concentration Calculation in ppb")
    print("This is calculated as: (moles of CH4 / total moles of atmosphere) * 1,000,000,000\n")
    
    # --- Final Output ---
    print("Final Equation:")
    print(f"ppb = (({final_ch4_mass_g:.2e} g / {mw_ch4} g/mol) / ({total_atmosphere_mass_g:e} g / {mw_air} g/mol)) * 1,000,000,000")
    print(f"ppb = ({final_ch4_moles:.2e} mol / {total_atmosphere_moles:.2e} mol) * 1,000,000,000")
    print(f"\nThe calculated increase in atmospheric methane concentration after 3 years is: {ppb_increase:.2f} ppb.")
    
    return ppb_increase

# Execute the calculation
final_ppb = calculate_methane_concentration_increase()
print(f"\n<<<80.25>>>")