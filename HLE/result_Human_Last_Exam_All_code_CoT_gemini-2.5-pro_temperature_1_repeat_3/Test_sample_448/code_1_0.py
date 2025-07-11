def calculate_rarest_noble_gas_percentage():
    """
    This script identifies the rarest noble gas and calculates its
    abundance as a percentage of terrestrial matter.

    The rarest noble gas is Radon (Rn). Its rarity is due to its
    radioactivity and short half-life. Its abundance is best estimated
    from its concentration in the Earth's crust.

    - Rarest Noble Gas: Radon (Rn)
    - Estimated abundance in Earth's crust (mass fraction): 4.0 x 10^-19
    """

    gas_name = "Radon (Rn)"
    # The abundance of Radon in the Earth's crust by mass fraction.
    abundance_mass_fraction = 4.0e-19

    # Convert the mass fraction to a percentage by multiplying by 100.
    percentage = abundance_mass_fraction * 100

    print(f"The rarest noble gas on Earth as of 2002 is {gas_name}.")
    print("Its abundance as a mass fraction in the Earth's crust is extremely low.")
    print("\nTo find the abundance as a percentage, we perform the following calculation:")
    
    # Using '{:e}' to format the numbers in scientific notation for clarity.
    print(f"{abundance_mass_fraction:e} (mass fraction) * 100 = {percentage:e}%")

calculate_rarest_noble_gas_percentage()