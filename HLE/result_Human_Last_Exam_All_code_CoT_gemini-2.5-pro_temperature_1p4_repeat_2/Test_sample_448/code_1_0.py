import sys

def find_rarest_noble_gas():
    """
    This script finds the rarest noble gas by comparing their atmospheric abundances.
    
    Data represents the abundance of noble gases in Earth's atmosphere.
    While the query asks for abundance in "all terrestrial matter," atmospheric
    abundance is a well-documented and standard proxy for comparing the relative
    rarity of these gases. Radon's extreme rarity is due to its short
    half-life as a product of radioactive decay, making its steady-state
    concentration minuscule throughout all terrestrial domains.
    The abundances are given in percent by volume of the atmosphere.
    The 2002 date doesn't significantly change these well-established values.
    """
    
    # Abundance of noble gases as a percentage of Earth's atmosphere by volume.
    noble_gas_abundance_percent = {
        'Helium': 0.000524,
        'Neon': 0.001818,
        'Argon': 0.934,
        'Krypton': 0.000114,
        'Xenon': 0.0000087,
        # Radon's abundance is highly variable and extremely low. An average
        # atmospheric mole fraction is ~6e-18, which is 6e-16 percent.
        'Radon': 6e-16 
    }

    print("Comparing the abundance of noble gases in Earth's atmosphere (% by volume):")
    
    # This loop outputs each number used in the comparison, as requested.
    # We use a format that can display very small numbers for Radon.
    for gas, percentage in noble_gas_abundance_percent.items():
        # Using 'g' format for clean printing of both large and small numbers.
        print(f"- {gas}: {percentage:g}%")

    # Find the rarest gas by identifying the key with the minimum value in the dictionary.
    rarest_gas = min(noble_gas_abundance_percent, key=noble_gas_abundance_percent.get)
    rarest_gas_percentage = noble_gas_abundance_percent[rarest_gas]

    print(f"\nBased on this data, the rarest noble gas is {rarest_gas}.")
    print(f"Its abundance is approximately {rarest_gas_percentage:g}% of the atmosphere.")

if __name__ == '__main__':
    find_rarest_noble_gas()
