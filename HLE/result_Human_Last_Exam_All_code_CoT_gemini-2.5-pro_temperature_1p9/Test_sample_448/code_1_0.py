import sys

def find_rarest_noble_gas():
    """
    Identifies the rarest noble gas on Earth based on its atmospheric abundance.
    The abundances are given as a percentage of all terrestrial matter,
    with atmospheric concentration being the primary contributor.
    """
    # Data represents abundance by volume percentage in Earth's atmosphere.
    # The abundance in "all terrestrial matter" is overwhelmingly dominated
    # by atmospheric concentration for these gases.
    noble_gas_abundance_percent = {
        'Argon (Ar)': 0.934,
        'Neon (Ne)': 0.001818,
        'Helium (He)': 0.000524,
        'Krypton (Kr)': 0.000114,
        'Xenon (Xe)': 0.0000087,
        'Radon (Rn)': 6e-18 # 0.000000000000000006%
    }

    # Find the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundance_percent, key=noble_gas_abundance_percent.get)
    rarest_gas_percentage = noble_gas_abundance_percent[rarest_gas_name]

    print("Comparing the abundance of noble gases as a percentage of terrestrial matter (primarily atmospheric):")
    for gas, percent in noble_gas_abundance_percent.items():
        # Use scientific notation for very small numbers for better readability
        print(f"- {gas}: {percent:.2e}%")

    print("\nThe rarest noble gas is the one with the lowest percentage.")
    sys.stdout.write("The rarest noble gas on Earth is Radon (Rn) with an approximate percentage of ")
    # The final print needs to show each number in the equation.
    # So we print the numbers of the percentage one by one.
    percentage_str = f"{rarest_gas_percentage:.1e}" # format to 6.0e-18
    parts = percentage_str.split('e')
    base = parts[0]
    exponent = parts[1]

    for char in base:
        sys.stdout.write(char)
    sys.stdout.write("e")
    for char in exponent:
        sys.stdout.write(char)
    print(" %.")


find_rarest_noble_gas()

print("<<<Radon>>>")