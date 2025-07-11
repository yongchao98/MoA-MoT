import math

def calculate_methane_increase():
    """
    Calculates the increase in atmospheric methane concentration after a 3-year period
    following the Nord Stream pipeline leak.
    """

    # --- Given Constants ---
    total_ch4_mass_metric_tons = 250000.0
    molar_mass_ch4 = 16.0  # g/mol
    mass_atmosphere_kg = 5.1e18
    molar_mass_air = 29.0  # g/mol, average molar mass of air

    # --- Step-by-Step Calculation ---

    # 1. Calculate initial moles of CH4 released.
    # Convert metric tons to grams (1 metric ton = 1e6 grams)
    total_ch4_mass_g = total_ch4_mass_metric_tons * 1e6
    initial_moles_ch4 = total_ch4_mass_g / molar_mass_ch4

    print("Step 1: Calculate initial moles of CH4 released.")
    print(f"Total CH4 Mass = {total_ch4_mass_metric_tons} metric tons")
    print(f"Initial Moles of CH4 = {total_ch4_mass_g:.2e} g / {molar_mass_ch4} g/mol = {initial_moles_ch4:.4e} moles")
    print("-" * 50)

    # 2. Calculate CH4 remaining after 3 years of mixing and oxidation.
    # Year 1: 80% in troposphere oxidized by 5%, 20% remains untouched initially.
    # The total fraction remaining is (0.8 * 0.95) + 0.2 = 0.76 + 0.2 = 0.96
    factor_y1 = (0.8 * (1 - 0.05)) + 0.2
    moles_after_y1 = initial_moles_ch4 * factor_y1

    # Year 2 & 3: Total reduced by 3% each year.
    factor_y2_y3 = (1 - 0.03)
    moles_after_y2 = moles_after_y1 * factor_y2_y3
    final_moles_ch4 = moles_after_y2 * factor_y2_y3

    print("Step 2: Calculate CH4 remaining after 3 years.")
    print(f"Moles after Year 1 reduction = {initial_moles_ch4:.4e} moles * {factor_y1:.2f} = {moles_after_y1:.4e} moles")
    print(f"Moles after Year 2 reduction = {moles_after_y1:.4e} moles * {factor_y2_y3:.2f} = {moles_after_y2:.4e} moles")
    print(f"Final Moles CH4 after Year 3 = {moles_after_y2:.4e} moles * {factor_y2_y3:.2f} = {final_moles_ch4:.4e} moles")
    print("-" * 50)

    # 3. Calculate total moles of air in the atmosphere.
    # Convert kg to g
    mass_atmosphere_g = mass_atmosphere_kg * 1e3
    total_moles_air = mass_atmosphere_g / molar_mass_air

    print("Step 3: Calculate total moles of air in the atmosphere.")
    print(f"Total Atmosphere Mass = {mass_atmosphere_kg:.2e} kg")
    print(f"Total Moles of Air = {mass_atmosphere_g:.2e} g / {molar_mass_air} g/mol = {total_moles_air:.4e} moles")
    print("-" * 50)

    # 4. Calculate the final increase in concentration in ppb.
    # ppb = (moles_substance / total_moles_air) * 1e9
    increase_in_ppb = (final_mles_ch4 / total_moles_air) * 1e9

    print("Step 4: Calculate the final increase in concentration (ppb).")
    print("Final Equation: (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
    print(f"Increase (ppb) = ({final_moles_ch4:.4e} / {total_moles_air:.4e}) * 1e9")
    print(f"Result: {increase_in_ppb:.4f} ppb")

    return increase_in_ppb

if __name__ == '__main__':
    final_answer = calculate_methane_increase()
    # The final answer is wrapped in <<<>>> as requested.
    print(f"<<<{final_answer:.4f}>>>")