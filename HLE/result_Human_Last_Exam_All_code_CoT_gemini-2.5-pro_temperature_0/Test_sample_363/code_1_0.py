def get_product_name():
    """
    This function returns the IUPAC name of the product from the described reaction.

    The reaction is a stereoselective aza-Ireland-Claisen rearrangement.
    The stereochemistry of the product is determined by the following:
    1. The starting material contains two chiral centers: (S) on the 1-phenylethyl group
       and (S) on the 5-methylcyclopentenyl group.
    2. The (S)-1-phenylethyl group acts as a chiral auxiliary, directing the formation
       of the new stereocenter on the propionamide chain to have an (R) configuration.
       This new center is at position 2 of the propionamide.
    3. The rearrangement creates a second new stereocenter on the cyclopentane ring.
       Its stereochemistry is determined relative to the other centers via a
       chair-like transition state, resulting in a (2S, 4S) configuration for the ring.
       The original (S) center at C5 becomes an (S) center at C4.
    4. The rearrangement also moves the double bond to an exocyclic position (a methylene group).
    """
    # The full IUPAC name of the final product amide.
    product_iupac_name = "(2R)-N-((S)-1-phenylethyl)-2-((2S,4S)-4-methyl-1-methylenecyclopentan-2-yl)propanamide"
    return product_iupac_name

# Print the final IUPAC name, including all numbers (locants).
final_product = get_product_name()
print("The IUPAC name of the product is:")
print(final_product)