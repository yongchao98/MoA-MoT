import sys

def find_rarest_noble_gas():
    """
    This script calculates the abundance of noble gases as a percentage of Earth's total mass
    to determine which one is the rarest.
    """
    # --- Constants ---
    # Mass of Earth in kg
    MASS_EARTH_KG = 5.972e24
    # Mass of Earth's atmosphere in kg
    MASS_ATMOSPHERE_KG = 5.15e18
    # Average molar mass of dry air in g/mol
    AVG_MOLAR_MASS_AIR_G_MOL = 28.97

    # Data for noble gases:
    # Abundance is in parts per million (ppm) by volume in the atmosphere.
    # Molar mass is in g/mol.
    # Note: Radon's concentration is extremely low and variable. 6e-14 ppm is a representative value
    # for the entire atmosphere, corresponding to a total mass of a few kilograms.
    noble_gases = [
        {'name': 'Helium (He)', 'abundance_ppm_v': 5.24, 'molar_mass': 4.0026},
        {'name': 'Neon (Ne)', 'abundance_ppm_v': 18.18, 'molar_mass': 20.180},
        {'name': 'Argon (Ar)', 'abundance_ppm_v': 9340, 'molar_mass': 39.948},
        {'name': 'Krypton (Kr)', 'abundance_ppm_v': 1.14, 'molar_mass': 83.798},
        {'name': 'Xenon (Xe)', 'abundance_ppm_v': 0.087, 'molar_mass': 131.29},
        {'name': 'Radon (Rn)', 'abundance_ppm_v': 6e-14, 'molar_mass': 222.0} # Extremely rare due to radioactivity
    ]

    rarest_gas = None
    min_percentage = float('inf')

    print("Calculating the percentage of each noble gas relative to total Earth mass...\n")

    for gas in noble_gases:
        # Step 1: Convert abundance from ppm (by volume) to a direct volume fraction
        volume_fraction = gas['abundance_ppm_v'] / 1_000_000

        # Step 2: Estimate mass fraction in the atmosphere
        # Mass fraction ≈ Volume fraction * (Molar Mass of Gas / Molar Mass of Air)
        mass_fraction_in_atmosphere = volume_fraction * (gas['molar_mass'] / AVG_MOLAR_MASS_AIR_G_MOL)

        # Step 3: Calculate the total mass of the gas in the atmosphere
        total_mass_of_gas_kg = mass_fraction_in_atmosphere * MASS_ATMOSPHERE_KG

        # Step 4: Calculate the percentage of the gas relative to Earth's total mass
        percentage_of_earth_mass = (total_mass_of_gas_kg / MASS_EARTH_KG) * 100

        if percentage_of_earth_mass < min_percentage:
            min_percentage = percentage_of_earth_mass
            rarest_gas = gas

    # Output the final result and show the calculation for the rarest gas
    if rarest_gas:
        print(f"The rarest noble gas on Earth is {rarest_gas['name']}.")
        print("\n--- Calculation Breakdown for the Rarest Gas ---")
        
        # Redo calculations for the final print statement
        rg_vol_frac = rarest_gas['abundance_ppm_v'] / 1_000_000
        rg_mass_frac_atm = rg_vol_frac * (rarest_gas['molar_mass'] / AVG_MOLAR_MASS_AIR_G_MOL)
        rg_mass_kg = rg_mass_frac_atm * MASS_ATMOSPHERE_KG
        rg_perc_earth = (rg_mass_kg / MASS_EARTH_KG) * 100

        print(f"1. Volume Fraction in Atmosphere: {rarest_gas['abundance_ppm_v']} / 1,000,000 = {rg_vol_frac:.2e}")
        print(f"2. Mass Fraction in Atmosphere: {rg_vol_frac:.2e} * ({rarest_gas['molar_mass']} g/mol / {AVG_MOLAR_MASS_AIR_G_MOL} g/mol) = {rg_mass_frac_atm:.2e}")
        print(f"3. Total Mass in Atmosphere (kg): {rg_mass_frac_atm:.2e} * {MASS_ATMOSPHERE_KG:.2e} kg = {rg_mass_kg:.2e} kg")
        print(f"4. Percentage of Earth's Mass: ({rg_mass_kg:.2e} kg / {MASS_EARTH_KG:.2e} kg) * 100 = {rg_perc_earth:.2e}%")
        
        print(f"\nThe percentage of {rarest_gas['name']} as a part of all terrestrial matter is approximately {min_percentage:.2e}%.")

if __name__ == '__main__':
    find_rarest_noble_gas()
