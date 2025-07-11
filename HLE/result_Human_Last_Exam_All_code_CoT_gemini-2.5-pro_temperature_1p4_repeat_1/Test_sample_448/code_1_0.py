def find_rarest_noble_gas():
    """
    This script identifies the rarest noble gas on Earth by comparing their
    abundances as a percentage of volume in the atmosphere. These values are
    a reliable proxy for overall rarity as the atmosphere is the primary
    reservoir for these elements. The abundances for 2002 are effectively
    the same as today's values for this comparison.
    """
    # Dictionary of noble gases and their abundance in Earth's atmosphere (in percent).
    # Radon's value is an estimate based on its typical atomic fraction.
    abundances_percent = {
        "Helium": 0.000524,
        "Neon": 0.001818,
        "Argon": 0.934,
        "Krypton": 0.000114,
        "Xenon": 0.0000087,
        "Radon": 6e-18 # (approx. 6.0 x 10^-18 %)
    }

    print("To find the rarest noble gas, we compare their percentages in Earth's atmosphere:")
    # The 'equation' is the comparison of all the values to find the minimum.
    # The following lines output each number involved in this comparison.
    for gas, percentage in abundances_percent.items():
        print(f"{gas}: {percentage:.2e} %")

    # Find the gas with the minimum percentage
    rarest_gas = min(abundances_percent, key=abundances_percent.get)
    min_percentage = abundances_percent[rarest_gas]

    print("\n---")
    print(f"The calculation shows that the gas with the lowest percentage is {rarest_gas}.")
    print(f"The final answer is that the rarest noble gas is {rarest_gas} with a percentage of approximately {min_percentage:.2e} %.")


find_rarest_noble_gas()