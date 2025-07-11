def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years.
    """
    # --- Step 1: Initial values and constants ---
    total_ch4_release_tons = 250000
    # Convert metric tons to grams (1 metric ton = 1,000,000 grams)
    initial_ch4_mass_g = total_ch4_release_tons * 1e6
    
    ch4_molar_mass_g_per_mol = 16.0
    # A standard approximation for the average molar mass of air.
    air_molar_mass_g_per_mol = 29.0
    total_atmosphere_mass_kg = 5.1e18
    total_atmosphere_mass_g = total_atmosphere_mass_kg * 1000

    # --- Step 2: Calculate mass remaining after Year 1 ---
    # 80% mixes into the troposphere and is reduced by 5%
    mass_in_tropo_initial = initial_ch4_mass_g * 0.80
    oxidation_loss_y1 = mass_in_tropo_initial * 0.05
    mass_in_tropo_after_ox = mass_in_tropo_initial - oxidation_loss_y1

    # 20% remains unmixed for year 1
    mass_unmixed = initial_ch4_mass_g * 0.20
    
    # Total mass at the end of Year 1
    total_mass_end_y1 = mass_in_tropo_after_ox + mass_unmixed
    
    # --- Step 3: Calculate mass remaining after Years 2 & 3 ---
    # This total mass is now subject to a 3% reduction per year for 2 years.
    # The reduction factor for two consecutive 3% reductions is (1 - 0.03)^2
    oxidation_factor_y2_y3 = (1 - 0.03)**2
    final_ch4_mass_g = total_mass_end_y1 * oxidation_factor_y2_y3

    # --- Step 4: Convert final mass to ppb concentration ---
    final_ch4_moles = final_ch4_mass_g / ch4_molar_mass_g_per_mol
    total_air_moles = total_atmosphere_mass_g / air_molar_mass_g_per_mol
    
    # ppb is the mole fraction (moles of CH4 / moles of air) multiplied by 1 billion
    ppb_increase = (final_ch4_moles / total_air_moles) * 1e9

    # --- Step 5: Output the results and the final equation ---
    print("Calculation Steps:")
    print(f"1. Initial CH4 Mass Emitted: {initial_ch4_mass_g:.2e} g")
    print(f"2. Mass Remaining After Year 1 (80% portion reduced by 5%): {total_mass_end_y1:.2e} g")
    print(f"3. Final Mass Remaining After 3 Years (reduced by 3% for 2 more years): {final_ch4_mass_g:.2e} g")
    print("-" * 30)
    print("Final Concentration Calculation:")
    print(f"Final Moles of CH4: {final_ch4_moles:.2e} mol")
    print(f"Total Moles of Air in Atmosphere: {total_air_moles:.2e} mol")
    print("-" * 30)
    
    print("Final Equation:")
    print(f"ppb_increase = (Final Moles CH4 / Total Moles Air) * 1e9")
    print(f"ppb_increase = ({final_ch4_moles:.4e} / {total_air_moles:.4e}) * 1e9")
    print(f"\nResulting Increase in CH4 Concentration: {ppb_increase:.2f} ppb")
    
    # The final answer in the required format
    print(f"\n<<<{ppb_increase:.3f}>>>")

calculate_methane_increase()