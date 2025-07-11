def solve_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration in ppb after 3 years
    following the Nord Stream leak, considering atmospheric mixing and oxidation.
    """
    # --- Step 1: Define Constants from the problem ---
    initial_ch4_release_metric_tons = 250000
    ch4_molar_mass_g_mol = 16.0
    air_molar_mass_g_mol = 29.0  # Standard approximation for molar mass of dry air
    atmosphere_mass_kg = 5.1 * 10**18
    troposphere_mixing_fraction = 0.80
    y1_tropospheric_oxidation_rate = 0.05
    y2_y3_oxidation_rate = 0.03
    simulation_duration_years = 3

    # --- Step 2: Calculate Initial Mass Distribution ---
    # Convert total release from metric tons to kilograms
    total_ch4_mass_kg = initial_ch4_release_metric_tons * 1000
    # Mass entering troposphere in year 1
    ch4_to_troposphere_total_kg = total_ch4_mass_kg * troposphere_mixing_fraction
    # Total mass entering stratosphere over 3 years
    ch4_to_stratosphere_total_kg = total_ch4_mass_kg * (1 - troposphere_mixing_fraction)
    # Assume linear mixing into the stratosphere over 3 years
    ch4_to_stratosphere_per_year_kg = ch4_to_stratosphere_total_kg / simulation_duration_years

    # --- Step 3: Simulate Atmospheric Changes Over 3 Years ---
    ch4_mass_in_troposphere_kg = 0.0
    ch4_mass_in_stratosphere_kg = 0.0

    # --- Year 1 ---
    # Add 80% of methane to the troposphere and apply 5% oxidation
    ch4_mass_in_troposphere_kg += ch4_to_troposphere_total_kg
    ch4_mass_in_troposphere_kg *= (1 - y1_tropospheric_oxidation_rate)
    # Add 1/3 of the stratospheric portion; no oxidation in year 1
    ch4_mass_in_stratosphere_kg += ch4_to_stratosphere_per_year_kg
    
    # --- Year 2 ---
    # Apply 3% oxidation to existing tropospheric methane
    ch4_mass_in_troposphere_kg *= (1 - y2_y3_oxidation_rate)
    # Add the next portion to the stratosphere, then apply 3% oxidation to the new total
    ch4_mass_in_stratosphere_kg += ch4_to_stratosphere_per_year_kg
    ch4_mass_in_stratosphere_kg *= (1 - y2_y3_oxidation_rate)

    # --- Year 3 ---
    # Apply 3% oxidation to existing tropospheric methane
    ch4_mass_in_troposphere_kg *= (1 - y2_y3_oxidation_rate)
    # Add the final portion to the stratosphere, then apply 3% oxidation to the new total
    ch4_mass_in_stratosphere_kg += ch4_to_stratosphere_per_year_kg
    ch4_mass_in_stratosphere_kg *= (1 - y2_y3_oxidation_rate)

    # --- Step 4: Sum Remaining Methane Mass ---
    total_remaining_ch4_mass_kg = ch4_mass_in_troposphere_kg + ch4_mass_in_stratosphere_kg

    # --- Step 5: Convert Mass to Moles ---
    # Convert molar masses from g/mol to kg/mol for consistency
    ch4_molar_mass_kg_mol = ch4_molar_mass_g_mol / 1000
    air_molar_mass_kg_mol = air_molar_mass_g_mol / 1000
    # Calculate moles of remaining methane
    moles_ch4_remaining = total_remaining_ch4_mass_kg / ch4_molar_mass_kg_mol
    # Calculate total moles of air in the atmosphere
    total_moles_air = atmosphere_mass_kg / air_molar_mass_kg_mol

    # --- Step 6: Calculate Final ppb Increase ---
    ppb_increase = (moles_ch4_remaining / total_moles_air) * 1e9

    # --- Print Results ---
    print("Final Calculation Steps:")
    print(f"1. The total remaining methane mass after 3 years is {total_remaining_ch4_mass_kg:.2f} kg.")
    print(f"2. Moles of remaining CH₄ = {moles_ch4_remaining:.2f} mol.")
    print(f"3. Total moles of air in the atmosphere = {total_moles_air:.2e} mol.")
    
    print("\nFinal Equation for ppb Increase:")
    print(f"ppb increase = (moles of remaining CH₄ / total moles of air) * 1,000,000,000")
    print(f"ppb increase = ({moles_ch4_remaining:.2f} / {total_moles_air:.2e}) * 1,000,000,000")
    print(f"\nThe calculated increase in atmospheric methane concentration after 3 years is {ppb_increase:.4f} ppb.")

    # Print the final answer in the required format
    print(f"\n<<<{ppb_increase:.4f}>>>")

solve_methane_concentration()