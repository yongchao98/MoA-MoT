import sys

def find_rarest_noble_gas():
    """
    This script identifies the rarest noble gas on Earth based on its
    atmospheric abundance and calculates its percentage.
    """
    # Noble gas abundances in Earth's atmosphere as a percentage by volume.
    # The data for 2002 is virtually identical to current standard values for this purpose.
    # Radon's abundance is an estimate (approx. 1 part in 10^21), as it varies and is extremely low.
    noble_gas_abundance_percent = {
        'Helium': 5.24e-4,
        'Neon': 1.818e-3,
        'Argon': 0.934,
        'Krypton': 1.14e-4,
        'Xenon': 8.7e-6,
        'Radon': 1e-19
    }

    # Find the gas with the minimum abundance.
    rarest_gas_name = min(noble_gas_abundance_percent, key=noble_gas_abundance_percent.get)
    rarest_gas_percentage = noble_gas_abundance_percent[rarest_gas_name]

    # The abundance of Radon is estimated at 1 part per 10^21.
    # The calculation for its percentage is (1 / 10^21) * 100.
    part = 1
    total_parts = "10^21" # Represented as a string for clear printing
    percentage_multiplier = 100

    print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
    print(f"Its approximate atmospheric abundance is {rarest_gas_percentage}%.")
    print(f"\nThis percentage is calculated from its estimated abundance of 1 part per 10^21:")
    # Using scientific notation for the result for readability
    print(f"Equation: ({part} / {total_parts}) * {percentage_multiplier} = {rarest_gas_percentage:.1e}%")

if __name__ == '__main__':
    find_rarest_noble_gas()