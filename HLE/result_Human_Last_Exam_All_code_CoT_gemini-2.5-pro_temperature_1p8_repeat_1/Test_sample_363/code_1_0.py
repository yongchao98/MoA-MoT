import sys

def get_product_name():
    """
    This function returns the IUPAC name of the product from the described reaction.

    The reaction is an Aza-Claisen ([3,3]-sigmatropic) rearrangement.
    1. LiHMDS creates a lithium enolate at the alpha-position of the propionamide.
    2. Heating causes the N-allyl group to rearrange to the alpha-carbon.
    3. The stereocenters in the starting material control the stereochemistry of the new chiral center created during the reaction.

    Starting Material: N-(((S)-5-methylcyclopent-1-en-1-yl)methyl)-N-((S)-1-phenylethyl)propionamide
    Predicted Product: A new amide where the allyl group has migrated and rearranged.
    """
    # The full IUPAC name for the major diastereomer formed.
    # Note the numbers required by the prompt are part of the name:
    # N-((S)-1-phenylethyl)-(2S)-2-((3S)-3-methyl-2-methylidenecyclopentan-1-yl)propanamide
    # Numbers present: 1, 2, 3, 2, 1
    product_name = "N-((S)-1-phenylethyl)-(2S)-2-((3S)-3-methyl-2-methylidenecyclopentan-1-yl)propanamide"
    return product_name

if __name__ == '__main__':
    final_product = get_product_name()
    print(final_product)