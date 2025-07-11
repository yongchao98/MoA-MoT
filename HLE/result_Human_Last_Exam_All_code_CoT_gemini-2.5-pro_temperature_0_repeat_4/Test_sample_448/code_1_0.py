def find_rarest_noble_gas():
    """
    Identifies the rarest noble gas on Earth by finding the minimum
    value among their atmospheric abundances. The year 2002 does not
    significantly change these established natural abundance values.
    """
    # A dictionary of noble gases and their abundance in Earth's atmosphere as a percentage.
    # "All terrestrial matter" is best represented by the atmosphere for noble gases.
    noble_gas_abundances = {
        'Helium (He)': 0.000524,
        'Neon (Ne)': 0.001818,
        'Argon (Ar)': 0.934,
        'Krypton (Kr)': 0.000114,
        'Xenon (Xe)': 0.0000087,
        # Radon (Rn) is radioactive, so its concentration is extremely low and variable.
        # The average is about 6x10^-18 mole percent, which is 6e-20 percent.
        'Radon (Rn)': 6e-20
    }

    print("To find the rarest noble gas, we must find the minimum value among their atmospheric percentages:")

    # Create a list of the percentage values to display in an "equation"
    abundance_values = list(noble_gas_abundances.values())

    # Find the minimum abundance value
    min_abundance = min(abundance_values)

    # Find the name of the gas corresponding to the minimum value
    rarest_gas_name = ""
    for gas, abundance in noble_gas_abundances.items():
        if abundance == min_abundance:
            rarest_gas_name = gas
            break

    # Format and print the comparison as a final equation
    # The numbers in the list correspond to He, Ne, Ar, Kr, Xe, and Rn respectively.
    equation_str = "min(" + ", ".join(f"{v:.7f}" if v > 1e-10 else f"{v:.1e}" for v in abundance_values) + f") = {min_abundance:.1e}"
    print(equation_str)

    print(f"\nThe minimum value in the equation, {min_abundance:.1e}, corresponds to {rarest_gas_name}.")
    print(f"Therefore, the rarest noble gas on Earth is {rarest_gas_name}, with an approximate percentage of {min_abundance:.1e}%.")

find_rarest_noble_gas()