import collections

def find_rarest_noble_gas():
    """
    This function identifies the rarest noble gas on Earth based on its
    atmospheric abundance, which is a proxy for its percentage in all terrestrial matter.
    The data is valid for 2002 as these concentrations are stable over long periods.
    """
    
    # Abundances of noble gases in Earth's atmosphere by percentage.
    # Sources: NASA, CRC Handbook of Chemistry and Physics.
    # Radon's abundance is exceptionally low due to its radioactivity and short half-life.
    noble_gas_percentages = {
        'Helium': 5.24e-4,      # 0.000524%
        'Neon': 1.818e-3,       # 0.001818%
        'Argon': 0.934,         # 0.934%
        'Krypton': 1.14e-4,     # 0.000114%
        'Xenon': 8.7e-6,        # 0.0000087%
        'Radon': 1e-19          # ~0.0000000000000000001%
    }

    # Find the gas with the minimum percentage
    rarest_gas_name = min(noble_gas_percentages, key=noble_gas_percentages.get)
    rarest_gas_percentage = noble_gas_percentages[rarest_gas_name]

    print("Comparison of Noble Gas Abundances (as percentage in atmosphere):")
    for gas, percentage in noble_gas_percentages.items():
        print(f"- {gas}: {percentage}%")
        
    print("\nResult:")
    # The final "equation" is the statement identifying the rarest gas and its value.
    print(f"The rarest noble gas is {rarest_gas_name} with an abundance of approximately {rarest_gas_percentage}%.")

if __name__ == '__main__':
    find_rarest_noble_gas()