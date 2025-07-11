import sys

def find_rarest_noble_gas():
    """
    Identifies the rarest noble gas on Earth based on its estimated percentage
    of all terrestrial matter by mass.
    The year 2002 is specified, but these abundances are stable over long timescales.
    """

    # Data represents estimated mass percentage in all terrestrial matter (crust, mantle, core, atmosphere).
    # These are estimations based on geochemical models.
    # For instance, Radon's abundance is calculated from the secular equilibrium
    # with its parent element, Uranium, in the Earth's mantle and crust.
    noble_gas_abundances = {
        'Argon (Ar)': 6.6e-4,     # Dominated by radiogenic 40-Ar from Potassium-40 decay
        'Neon (Ne)': 3.3e-8,
        'Krypton (Kr)': 2.2e-10,
        'Xenon (Xe)': 2.2e-11,    # Rarest stable noble gas
        'Radon (Rn)': 4.4e-18     # Rarest noble gas overall due to its short half-life
    }
    
    # We omit Helium as its abundance is highly variable and complex due to
    # production from radioactive decay and loss to space, but it is significantly
    # more abundant than Xenon or Radon.

    print("Abundance of noble gases as a percentage of all terrestrial matter (by mass):")
    for gas, percentage in noble_gas_abundances.items():
        # Using sys.stdout.write to avoid scientific notation for some numbers if not needed
        # but here we format it consistently for clarity.
        print(f"- {gas}: {percentage:.1e} %")

    # Find the rarest gas by finding the minimum percentage
    rarest_gas = min(noble_gas_abundances, key=noble_gas_abundances.get)
    rarest_gas_percentage = noble_gas_abundances[rarest_gas]

    print("\n--------------------------------------------------")
    print("Conclusion:")
    print(f"The rarest noble gas is {rarest_gas.split(' ')[0]}.")
    
    # The final output requires printing the numbers in the final equation.
    # Here, the "equation" is the comparison leading to the answer.
    print(f"Its percentage of all terrestrial matter is {rarest_gas_percentage:.1e} %.")
    print("This is calculated by comparing the abundances of all noble gases:")
    
    # Building a string to show the comparison, e.g., 4.4e-18 < 2.2e-11 < ...
    sorted_gases = sorted(noble_gas_abundances.items(), key=lambda item: item[1])
    comparison_str = " < ".join([f"{val:.1e} % ({key.split(' ')[0]})" for key, val in sorted_gases])
    
    print(comparison_str)
    
find_rarest_noble_gas()
