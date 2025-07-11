def get_product_iupac_name():
    """
    This function returns the IUPAC name of the product from the given reaction.
    The reaction is a domino sequence: Oxy-Cope rearrangement, tautomerization,
    and a transannular Michael addition.
    """
    # Name components of the final product
    stereochemistry = "(E)"
    substituent_1_position = 5
    substituent_1_name = "methoxy"
    substituent_2_position = 6
    substituent_2_name = "methyl"
    parent_ring_system = "bicyclo[4.4.0]dec"
    unsaturation_position = 1
    unsaturation_type = "en"
    ketone_position = 4
    ketone_suffix = "one"

    # Constructing the IUPAC name string by string
    name = (
        f"{stereochemistry}-"
        f"{substituent_1_position}-{substituent_1_name}-"
        f"{substituent_2_position}-{substituent_2_name}"
        f"{parent_ring_system}-"
        f"{unsaturation_position}-{unsaturation_type}-"
        f"{ketone_position}-{ketone_suffix}"
    )
    
    # Print the equation part by part for clarity
    print(f"Final IUPAC Name: {stereochemistry} - {substituent_1_position}-methoxy - {substituent_2_position}-methyl - bicyclo[4.4.0]dec - {unsaturation_position}-en - {ketone_position}-one")
    
    # Print the final combined name
    print(f"The full IUPAC name is: {name}")


if __name__ == "__main__":
    get_product_iupac_name()