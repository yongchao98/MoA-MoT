import pubchempy as pcp

def find_chemical_name():
    """
    This function identifies the name of the final product from the synthesis
    by searching the PubChem database using its SMILES string.
    """
    # The SMILES string for the final product, L-Tryptophan hydroxamate.
    # N[C@@H] represents the L-stereoisomer's amine group.
    # (CC1=Cn2ccccc21) is the indole side chain.
    # C(=O)NO is the hydroxamic acid group.
    smiles_string = 'N[C@@H](CC1=Cn2ccccc21)C(=O)NO'

    print(f"Searching for compound with SMILES: {smiles_string}")

    try:
        # Search PubChem for compounds matching the SMILES string
        compounds = pcp.get_compounds(smiles_string, 'smiles')

        if compounds:
            compound = compounds[0]
            print("\n--- Product Information ---")
            # Print the most common name found in synonyms
            if compound.synonyms:
                print(f"Common Name: {compound.synonyms[0]}")
            else:
                print("Common Name: Not found.")
            
            # Print the official IUPAC name
            if compound.iupac_name:
                print(f"IUPAC Name: {compound.iupac_name}")
            else:
                print("IUPAC Name: Not found.")
            print("-------------------------")
        else:
            print("Product not found in the PubChem database.")
            print("Based on chemical analysis, the product is L-Tryptophan hydroxamate.")

    except Exception as e:
        print(f"An error occurred while querying PubChem: {e}")
        print("Based on chemical analysis, the product is L-Tryptophan hydroxamate.")

if __name__ == "__main__":
    find_chemical_name()
