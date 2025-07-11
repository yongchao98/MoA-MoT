def find_rarest_noble_gas():
    """
    This function identifies the rarest stable noble gas in Earth's atmosphere
    and calculates its abundance as a percentage.
    """
    # Abundance of stable noble gases in Earth's atmosphere in parts per million (ppm)
    # Data is from 2002 or generally accepted values around that time.
    noble_gas_abundances_ppm = {
        'Helium': 5.24,
        'Neon': 18.18,
        'Argon': 9340,
        'Krypton': 1.14,
        'Xenon': 0.087
    }

    # Find the rarest gas by finding the minimum value in the dictionary
    rarest_gas_name = min(noble_gas_abundances_ppm, key=noble_gas_abundances_ppm.get)
    abundance_in_ppm = noble_gas_abundances_ppm[rarest_gas_name]

    # Convert abundance from ppm to percentage (1 ppm = 0.0001%)
    abundance_in_percent = abundance_in_ppm / 10000

    # Print the result and the calculation
    print(f"The rarest stable noble gas in Earth's atmosphere is {rarest_gas_name}.")
    print("The calculation for its percentage of all terrestrial matter (approximated by atmospheric abundance) is:")
    # Using f-string to format the output to show the full precision of the float
    print(f"{abundance_in_ppm} ppm / 10000 = {abundance_in_percent:.7f}%")

if __name__ == '__main__':
    find_rarest_noble_gas()