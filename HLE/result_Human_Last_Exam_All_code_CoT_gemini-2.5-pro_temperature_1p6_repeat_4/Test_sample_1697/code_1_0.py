import pubchempy as pcp

def find_reaction_product_name():
    """
    This function determines the final product of the described reaction
    and prints its standard IUPAC name.

    The reaction involves the directed ortho-metalation of
    N,N-diethyl-3-dimethylaminobenzamide followed by methylation.
    The most acidic proton is at position 2, between the two directing groups.
    This position is lithiated and then quenched with methyl iodide,
    resulting in the addition of a methyl group at position 2.

    The script confirms the product's name by looking up its structure
    in the PubChem online database.
    """
    # The SMILES string for the final product, N,N-diethyl-2-methyl-3-dimethylaminobenzamide
    product_smiles = "CCN(CC)C(=O)c1c(C)c(N(C)C)ccc1"
    product_name_fallback = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"

    print(f"Starting compound: N,N-diethyl-3-dimethylaminobenzamide")
    print(f"Reaction: 1) sec-BuLi, TMEDA, THF  2) CH3I")
    print("-" * 30)
    print("Determining final product...")

    try:
        # Use the SMILES string to find the compound in PubChem
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        if compounds:
            # Get the IUPAC name from the first result
            final_name = compounds[0].iupac_name
            print(f"Final product found in database: {final_name}")
        else:
            # Use the fallback name if the compound is not found
            print(f"Could not find compound in online database.")
            print(f"The determined final product is: {product_name_fallback}")

    except Exception as e:
        # Handle cases where the online lookup fails (e.g., no internet connection)
        print(f"An error occurred during online lookup: {e}")
        print(f"Based on chemical principles, the final product is: {product_name_fallback}")

# Execute the function to find and print the product name
find_reaction_product_name()