def find_rarest_noble_gas():
    """
    Identifies the rarest noble gas based on its atmospheric abundance and calculates its percentage.
    The year 2002 does not significantly change these natural abundances.
    """
    # Data for noble gas abundance in Earth's atmosphere (mole fraction).
    # Atmospheric abundance is a good proxy for overall terrestrial abundance.
    noble_gas_abundance = {
        'Helium': 5.24e-6,
        'Neon': 1.818e-5,
        'Argon': 9.34e-3,
        'Krypton': 1.14e-6,
        'Xenon': 8.7e-8,
        'Radon': 6e-18  # Source: Various chemistry resources, represents an average atmospheric mole fraction.
    }

    # Find the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundance, key=noble_gas_abundance.get)
    min_abundance_fraction = noble_gas_abundance[rarest_gas_name]

    # Convert the fraction to a percentage
    percentage = min_abundance_fraction * 100

    print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
    print("Its abundance as a percentage is calculated from its mole fraction in the atmosphere.")
    print("Final Equation:")
    # The format string "{:.1e}" displays the number in scientific notation with one decimal place.
    print(f"{min_abundance_fraction:.1e} (mole fraction) * 100 = {percentage:.1e}%")

if __name__ == '__main__':
    find_rarest_noble_gas()