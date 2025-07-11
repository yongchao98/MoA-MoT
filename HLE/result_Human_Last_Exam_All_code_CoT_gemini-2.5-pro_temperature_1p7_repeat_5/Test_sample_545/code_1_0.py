def solve_chemistry_problem():
    """
    This function determines the IUPAC name of the major product and prints the result.
    The logic follows a two-step reaction sequence:
    1. Thermal sulfoxide syn-elimination.
    2. Claisen rearrangement of the resulting allyl vinyl ether.
    """

    # Step 1: Define the numbers (locants) in the final IUPAC name.
    methyl_locant_1 = 5
    methyl_locant_2 = 5
    ene_locant = 4

    # Step 2: Define the parts of the name.
    prefix = "dimethyl"
    parent_chain = "pent"
    suffix_ene = "en"
    suffix_aldehyde = "al"

    # Step 3: Construct the final IUPAC name.
    final_name = (f"{methyl_locant_1},{methyl_locant_2}-{prefix}"
                  f"{parent_chain}-{ene_locant}-{suffix_ene}{suffix_aldehyde}")

    # Step 4: Print the final answer, showing the numbers used.
    print("The IUPAC name of the major product is constructed as follows:")
    print(f"Methyl group locants: {methyl_locant_1}, {methyl_locant_2}")
    print(f"Alkene group locant: {ene_locant}")
    print("\nFinal IUPAC Name:")
    print(final_name)

solve_chemistry_problem()