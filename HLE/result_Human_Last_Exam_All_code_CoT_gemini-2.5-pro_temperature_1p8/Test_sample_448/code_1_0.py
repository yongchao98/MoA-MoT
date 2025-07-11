import sys

def find_rarest_noble_gas_percentage():
    """
    This script identifies the rarest noble gas and calculates its percentage
    of all terrestrial matter based on data from 2002.

    The rarest noble gas is Radon (Rn). Its rarity is due to its radioactivity
    and short half-life; it decays soon after being formed.

    Data is sourced from the CRC Handbook of Chemistry and Physics, 83rd Edition (2002-2003),
    which lists Radon's abundance in the Earth's crust. This is a standard proxy
    for 'terrestrial matter'.
    """

    # Rarest noble gas name
    gas_name = "Radon"

    # Abundance in Earth's crust in parts per million (ppm) by mass.
    # Source value is 4 x 10^-13 ppm.
    abundance_ppm = 4e-13

    # To convert parts per million (ppm) to a percentage, we divide by 10,000.
    # 1% = 10,000 ppm
    conversion_factor = 10000
    abundance_percentage = abundance_ppm / conversion_factor

    # To fulfill the request to output each number in the final equation,
    # we will represent the final percentage in the format: base x 10^exponent.
    # abundance_percentage = 4e-17
    base = 4.0
    exponent = -17

    print(f"The rarest noble gas on Earth as of 2002 is {gas_name}.")
    print(f"Its abundance as a percentage of terrestrial matter (crustal) is approximately:")
    print(f"{base} x 10^{exponent}%")

# Execute the function to print the result
find_rarest_noble_gas_percentage()