def find_rarest_noble_gas():
    """
    This script finds the rarest noble gas on Earth based on its
    atmospheric abundance and calculates its percentage.
    The abundances are given for 2002, but these values are stable over time.
    """
    # Abundance of noble gases in Earth's atmosphere in parts per million (ppm) by volume.
    # Data is a well-established scientific consensus.
    # Radon's abundance is extremely low (approx. 1 part in 10^21) due to its radioactivity.
    # 1 part in 10^21 is 1e-15 ppm.
    noble_gas_abundances_ppm = {
        'Helium': 5.24,
        'Neon': 18.18,
        'Argon': 9340,
        'Krypton': 1.14,
        'Xenon': 0.087,
        'Radon': 1e-15
    }

    # Find the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundances_ppm, key=noble_gas_abundances_ppm.get)
    rarest_gas_ppm = noble_gas_abundances_ppm[rarest_gas_name]

    # Convert the abundance from ppm to a percentage
    # Conversion factor: 1 ppm = 0.0001 %
    conversion_factor = 10000
    rarest_gas_percentage = rarest_gas_ppm / conversion_factor

    # Print the results
    print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
    print(f"The equation to find its percentage from its abundance in parts per million (ppm) is:")
    # Using scientific notation for better readability of the small numbers
    print(f"{rarest_gas_ppm:.1e} ppm / {conversion_factor} = {rarest_gas_percentage:.1e} %")
    print(f"As a percentage of all terrestrial matter (approximated by atmospheric abundance), this is {rarest_gas_percentage:.21f}%.")

find_rarest_noble_gas()