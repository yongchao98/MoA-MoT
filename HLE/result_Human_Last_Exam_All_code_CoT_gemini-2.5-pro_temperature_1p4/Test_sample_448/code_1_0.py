def find_rarest_noble_gas():
    """
    Identifies the rarest noble gas on Earth based on its atmospheric abundance.
    The data for 2002 is effectively the same as current data for these gases.
    """
    # A dictionary stores the noble gases and their approximate abundance
    # in Earth's atmosphere as a percentage by volume.
    noble_gas_abundances = {
        "Argon (Ar)": 0.934,
        "Neon (Ne)": 0.001818,
        "Helium (He)": 0.000524,
        "Krypton (Kr)": 0.000114,
        "Xenon (Xe)": 8.7e-6,  # 0.0000087%
        "Radon (Rn)": 1e-19    # Radon's abundance is extremely low, about 1 part per 10^21, which is 1e-19 percent.
    }

    # Find the gas with the minimum abundance from the dictionary
    rarest_gas_name = min(noble_gas_abundances, key=noble_gas_abundances.get)
    rarest_gas_percentage = noble_gas_abundances[rarest_gas_name]

    # Print the final result including the name and percentage value.
    print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
    print("Its approximate abundance as a percentage of the atmosphere is:")
    # Formatting the output to show the number in the final statement
    print(f"Abundance of {rarest_gas_name} = {rarest_gas_percentage:.1e}%")

if __name__ == '__main__':
    find_rarest_noble_gas()