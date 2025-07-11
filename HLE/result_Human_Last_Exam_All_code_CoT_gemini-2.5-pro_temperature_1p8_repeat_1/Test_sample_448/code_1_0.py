import sys

def find_rarest_noble_gas():
    """
    This script determines the rarest noble gas on Earth by comparing their
    abundances in the atmosphere.
    """
    # The abundance of noble gases in Earth's atmosphere as a percentage by volume.
    # These values are fundamentally constant and were the same in 2002.
    # Source: CRC Handbook of Chemistry and Physics & other scientific sources.
    noble_gas_abundances = {
        'Helium': 0.00052,
        'Neon': 0.0018,
        'Argon': 0.934,
        'Krypton': 0.000114,
        'Xenon': 0.0000087,
        # Radon is a radioactive gas with a highly variable and extremely low concentration.
        # Its atmospheric percentage is on the order of 10^-19 %, making it the rarest by far.
        'Radon': 1e-19
    }

    print("Step 1: Define the percentage of each noble gas in Earth's atmosphere.")
    for gas, percentage in noble_gas_abundances.items():
        # Using scientific notation for very small numbers for better readability
        print(f"{gas}: {percentage:.{1 if '%' in str(percentage) else '10g'}}%")

    # Find the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundances, key=noble_gas_abundances.get)
    min_abundance = noble_gas_abundances[rarest_gas_name]
    
    # Building the "equation" string as requested by showing the numbers to be compared.
    abundance_list = list(noble_gas_abundances.values())

    print("\nStep 2: Find the minimum value among these percentages.")
    print("The values being compared are:")
    # sys.stdout.write is used to avoid newlines for the "equation"
    sys.stdout.write("min(")
    for i, val in enumerate(abundance_list):
        sys.stdout.write(f"{val}")
        if i < len(abundance_list) - 1:
            sys.stdout.write(", ")
    sys.stdout.write(")\n")

    print(f"\nStep 3: The result of this comparison is {min_abundance}%.")
    print(f"\nThis corresponds to the noble gas: {rarest_gas_name}.")


if __name__ == '__main__':
    find_rarest_noble_gas()