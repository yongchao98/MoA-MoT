import sys

def find_rarest_noble_gas():
    """
    This script determines the rarest noble gas on Earth based on its
    atmospheric abundance.
    """
    # Atmospheric abundance of noble gases in parts per million (ppm).
    # Radon's abundance is highly variable and trace, but it is many orders
    # of magnitude lower than any other noble gas. We use a very small number
    # (1e-19, or 0.0000000000000000001 ppm) to represent its extreme rarity for calculation.
    # The abundance of Oganesson is zero as it is a synthetic element.
    noble_gas_abundances_ppm = {
        'Helium (He)': 5.24,
        'Neon (Ne)': 18.18,
        'Argon (Ar)': 9340,
        'Krypton (Kr)': 1.14,
        'Xenon (Xe)': 0.087,
        'Radon (Rn)': 1e-19,
    }

    print("Comparing the atmospheric abundance of noble gases in parts per million (ppm):")
    # This loop outputs each number used in the comparison.
    for gas, ppm in noble_gas_abundances_ppm.items():
        # Use scientific notation for very small or large numbers for readability
        if ppm < 1e-3 or ppm > 1e3:
            print(f"- {gas}: {ppm:e} ppm")
        else:
            print(f"- {gas}: {ppm} ppm")

    # Find the gas with the minimum abundance
    rarest_gas_name = min(noble_gas_abundances_ppm, key=noble_gas_abundances_ppm.get)
    rarest_gas_ppm = noble_gas_abundances_ppm[rarest_gas_name]

    # Convert ppm to percentage for the final answer. 1 ppm = 0.0001%
    rarest_gas_percentage = rarest_gas_ppm / 10000

    print("\nBased on this data, the rarest noble gas is Radon (Rn).")
    print("\nTo express its abundance as a percentage of the atmosphere, we perform the following calculation:")
    
    # Final equation showing the conversion from ppm to percentage.
    # We use '.{}' format to avoid forcing scientific notation on a potentially small number.
    print(f"Abundance Percentage = {rarest_gas_ppm} ppm / 10,000 = {rarest_gas_percentage}%")

if __name__ == '__main__':
    find_rarest_noble_gas()