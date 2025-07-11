def generate_iupac_name():
    """
    Constructs and prints the IUPAC name for the final reaction product.
    The final name is Ethyl 2,5-dihydrothiophene-3-carboxylate.
    This script will print the name and the numbers used to form it.
    """

    # Parts of the IUPAC name
    prefix = "Ethyl"
    locant_a = "2"
    locant_b = "5"
    parent_ring = "dihydrothiophene"
    locant_c = "3"
    suffix = "carboxylate"

    # Assemble the final name, printing each number as part of the "equation"
    print("The final name is constructed as follows:")
    final_name_eq = f"{prefix} {locant_a},{locant_b}-{parent_ring}-{locant_c}-{suffix}"

    print("Name construction: " + prefix + " " + locant_a + "," + locant_b + "-" + parent_ring + "-" + locant_c + "-" + suffix)
    print("\nThe complete IUPAC name is:")
    print(final_name_eq)

generate_iupac_name()