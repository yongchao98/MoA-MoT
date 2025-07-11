def find_rarest_noble_gas():
    """
    This function identifies the rarest noble gas on Earth based on its
    atmospheric abundance and calculates its percentage.
    """
    # The abundances of noble gases in Earth's atmosphere are a good proxy for their overall rarity.
    # The date 2002 doesn't matter, as these values are stable over human history.
    # Abundances are given in parts per million (ppm) by volume.
    # Note: Radon's concentration is extremely low and variable due to its short half-life (3.8 days).
    # Its concentration is often cited on the order of 1 part in 10^21, which corresponds to 1e-15 ppm.
    noble_gas_abundance_ppm = {
        'Helium': 5.24,
        'Neon': 18.18,
        'Argon': 9340,
        'Krypton': 1.14,
        'Xenon': 0.087,
        'Radon': 1e-15,  # An approximation of Radon's minuscule concentration
    }

    # Find the name of the gas with the lowest abundance
    rarest_gas_name = min(noble_gas_abundance_ppm, key=noble_gas_abundance_ppm.get)

    # Get the ppm value for the rarest gas
    rarest_gas_ppm = noble_gas_abundance_ppm[rarest_gas_name]

    # Convert the ppm value to a percentage
    # Percentage = (parts_per_million / 1,000,000) * 100
    percentage = (rarest_gas_ppm / 1_000_000) * 100

    print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
    print("Based on its atmospheric abundance, the value as a percentage is:")
    
    # We print the final result as an equation using the calculated numbers.
    # We use scientific notation for such a small number.
    ppm_val_str = f"{rarest_gas_ppm:.1e}"
    percentage_val_str = f"{percentage:.1e}"

    print(f"{rarest_gas_name} percentage = ({ppm_val_str} ppm / 1000000) * 100 = {percentage_val_str} %")

find_rarest_noble_gas()