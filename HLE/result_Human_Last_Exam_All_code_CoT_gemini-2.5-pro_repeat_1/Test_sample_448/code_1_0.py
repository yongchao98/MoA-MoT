import sys

def find_rarest_noble_gas():
    """
    This script identifies the rarest noble gas on Earth based on its
    abundance in terrestrial matter and calculates its percentage.
    """
    # Data on the abundance of noble gases in the Earth's crust by mass fraction.
    # Crustal abundance is a good proxy for comparing the rarity of these elements
    # in "all terrestrial matter". The values for stable elements do not change
    # significantly over decades, so 2002 data is effectively the same as current data.
    # Source: CRC Handbook of Chemistry and Physics, and other geological sources.
    noble_gas_abundances = {
        'Helium': 8e-9,      # 8 parts per billion
        'Neon': 7e-11,       # 70 parts per trillion
        'Argon': 3.5e-6,     # 3.5 parts per million
        'Krypton': 1e-10,    # 100 parts per trillion
        'Xenon': 3e-11,      # 30 parts per trillion
        'Radon': 4e-19       # 0.4 parts per quintillion (from radioactive decay)
    }

    # Find the rarest gas by identifying the minimum abundance value.
    # The 'min' function with a 'key' argument finds the dictionary key
    # associated with the minimum value.
    rarest_gas_name = min(noble_gas_abundances, key=noble_gas_abundances.get)
    min_abundance_fraction = noble_gas_abundances[rarest_gas_name]

    # Convert the mass fraction to a percentage by multiplying by 100.
    abundance_percentage = min_abundance_fraction * 100

    # Print the final result in a clear format.
    print(f"The rarest noble gas on Earth as of 2002 (and today) is {rarest_gas_name}.")
    print("Its abundance as a percentage of all terrestrial matter is calculated as follows:")

    # Print the calculation as requested, showing the numbers involved.
    # Using scientific notation for clarity due to the extremely small value.
    print(f"Mass Fraction: {min_abundance_fraction}")
    print(f"Percentage = {min_abundance_fraction} * 100")
    print(f"Result: {abundance_percentage:.1e} %")

if __name__ == '__main__':
    find_rarest_noble_gas()