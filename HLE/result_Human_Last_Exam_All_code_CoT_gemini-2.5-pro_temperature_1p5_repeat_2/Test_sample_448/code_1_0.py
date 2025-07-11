import sys

def find_rarest_noble_gas():
    """
    This function identifies the rarest noble gas on Earth based on its
    abundance in the atmosphere, which serves as a proxy for its overall rarity.
    """

    # Data representing the abundance of noble gases in Earth's atmosphere as a percentage.
    # The values from 2002 are not materially different from these standard figures.
    noble_gas_abundance_percent = {
        'Argon': 0.934,
        'Neon': 0.001818,
        'Helium': 0.000524,
        'Krypton': 0.000114,
        'Xenon': 0.0000087,
        # Radon is radioactive and decays quickly, so it exists only in trace amounts.
        # Its average concentration is orders of magnitude smaller than Xenon's.
        # A value of 1e-19 percent is a representative, extremely small figure.
        'Radon': 1e-19
    }

    # Find the gas with the minimum abundance
    # The min() function's key argument is used to specify a function to be called on each
    # list element prior to making comparisons. Here, it gets the value from the dictionary.
    rarest_gas_name = min(noble_gas_abundance_percent, key=noble_gas_abundance_percent.get)

    # Retrieve the abundance value for the rarest gas
    rarest_gas_percentage = noble_gas_abundance_percent[rarest_gas_name]

    # Print the final result, showing the "equation" of which gas and percentage is the answer.
    # The f-string formats the very small number into scientific notation for readability.
    print(f"To find the rarest noble gas, we compare their abundances:")
    for gas, percent in noble_gas_abundance_percent.items():
        print(f"- {gas}: {percent}%")
    
    print("\nAfter comparison, the conclusion is:")
    print(f"The rarest noble gas = {rarest_gas_name}")
    print(f"Its percentage = {rarest_gas_percentage:.1e}%")


find_rarest_noble_gas()