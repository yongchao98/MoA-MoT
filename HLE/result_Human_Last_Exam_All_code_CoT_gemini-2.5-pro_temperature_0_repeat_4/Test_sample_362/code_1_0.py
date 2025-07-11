def solve_wittig_reaction():
    """
    This script determines the product of a Wittig reaction by identifying the
    reactants, breaking them down into their reacting fragments, and combining
    the fragments to form the final product.
    """

    # 1. Define the reactants
    aldehyde_name = "pivalaldehyde"
    aldehyde_structure = "(CH3)3C-CHO"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_structure = "(2-Cl-C6H4)-CH2-CH=PPh3"

    # 2. Identify the fragments that will be joined
    # The aldehyde contributes the part attached to the carbonyl group
    aldehyde_fragment = "(CH3)3C-CH"
    # The ylide contributes the carbon part of the C=P bond
    ylide_fragment = "CH-CH2-(2-Cl-C6H4)"

    # 3. Construct the product and byproduct structures
    # The C=O and C=P are replaced by a new C=C bond
    product_structure = f"{aldehyde_fragment}={ylide_fragment}"
    byproduct_structure = "O=PPh3"

    # 4. Determine the IUPAC name of the product
    # The main chain is a pentene (5 carbons).
    # Numbering from the phenyl side gives the double bond the lowest number (2).
    # 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene
    product_iupac_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # 5. Print the full reaction equation and explain the numbers in the name
    print("The Wittig reaction is:")
    print(f"{aldehyde_structure}  +  {ylide_structure}  --->  {product_structure}  +  {byproduct_structure}\n")
    print(f"The main organic product is named: {product_iupac_name}\n")
    print("Explanation of the numbers in the IUPAC name:")
    print("The final equation for the product name is built from its components:")
    print("Position '1': The (2-chlorophenyl) group is on the first carbon of the main chain.")
    print("Positions '4,4': Two methyl groups are on the fourth carbon of the main chain.")
    print("Position '2': The double bond ('-ene') starts at the second carbon of the main chain.")

solve_wittig_reaction()
<<<1-(2-chlorophenyl)-4,4-dimethylpent-2-ene>>>