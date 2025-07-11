import sys

def find_rarest_noble_gas():
    """
    This function identifies the rarest stable noble gas from a predefined list
    of their atmospheric abundances, and calculates its percentage.
    """
    # Step 1: Define the abundance of stable noble gases in Earth's atmosphere
    # in parts per million (ppm) by volume. These concentrations are stable over
    # decades, so the 2002 values are the same as today's accepted values.
    # Note: Radon is the absolute rarest noble gas, but it is radioactive with a
    # highly variable concentration. Xenon is the rarest stable noble gas.
    noble_gas_abundances_ppm = {
        'Helium': 5.24,
        'Neon': 18.18,
        'Argon': 9340,
        'Krypton': 1.14,
        'Xenon': 0.087
    }

    # Step 2: Find the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundances_ppm, key=noble_gas_abundances_ppm.get)
    rarest_gas_ppm = noble_gas_abundances_ppm[rarest_gas_name]

    # Step 3: Calculate the percentage from ppm
    # The conversion equation is: percentage = ppm / 10,000
    rarest_gas_percentage = rarest_gas_ppm / 10000

    # Step 4: Print the results, including the numbers in the final equation
    print(f"The rarest stable noble gas on Earth is {rarest_gas_name}.")
    print(f"Its atmospheric abundance is approximately {rarest_gas_ppm} parts per million (ppm).")
    print("\nThe equation to convert its abundance to a percentage is:")
    # Using 'g' format specifier to handle scientific notation for small numbers gracefully.
    print(f"Percentage = Abundance in ppm / 10000")
    print(f"Percentage = {rarest_gas_ppm} / 10000 = {rarest_gas_percentage:.9f}%")


if __name__ == '__main__':
    find_rarest_noble_gas()