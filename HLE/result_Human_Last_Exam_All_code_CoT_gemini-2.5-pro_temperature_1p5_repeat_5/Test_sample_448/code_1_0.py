def find_rarest_noble_gas():
    """
    This script identifies the rarest noble gas on Earth by comparing their
    atmospheric abundances, which serve as a reliable proxy for their
    concentration in all terrestrial matter. The data is valid for 2002.
    """

    # A dictionary mapping each noble gas to its abundance in percent (%) by volume
    # in Earth's atmosphere. Radon's concentration is an average due to its
    # radioactivity and high variability.
    # Sources: CRC Handbook of Chemistry and Physics, National Center for Atmospheric Research.
    noble_gas_abundance_percent = {
        'Argon': 9.34e-1,    # 0.934%
        'Neon': 1.818e-3,    # 0.001818%
        'Helium': 5.24e-4,     # 0.000524%
        'Krypton': 1.14e-4,    # 0.000114%
        'Xenon': 8.7e-6,     # 0.0000087%
        'Radon': 6.0e-18     # 0.000000000000000006%
    }

    # Find the name of the gas with the smallest percentage value
    rarest_gas_name = min(noble_gas_abundance_percent, key=noble_gas_abundance_percent.get)
    rarest_gas_percentage = noble_gas_abundance_percent[rarest_gas_name]

    # Output the final answer in a clear format, showing the name and the numeric value.
    print("The rarest noble gas on Earth and its percentage abundance is:")
    print(f"{rarest_gas_name} = {rarest_gas_percentage}%")

find_rarest_noble_gas()