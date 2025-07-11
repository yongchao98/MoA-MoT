def get_compound_4_details():
    """
    This function provides the identity of Compound 4, the final product of the reaction sequence.
    """
    # Names for the final compound
    common_name = "2,2'-dicarboxybenzophenone"
    iupac_name = "2,2'-Carbonylbis(benzoic acid)"

    # Chemical Formula
    molecular_formula = "C15H10O5"

    # A standard chemical structure representation (SMILES)
    smiles_string = "O=C(O)c1ccccc1C(=O)c2ccccc2C(=O)O"

    # Print the identifying information for Compound 4
    print("--- Identity of Final Product (Compound 4) ---")
    print(f"Common Name: {common_name}")
    print(f"IUPAC Name: {iupac_name}")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"SMILES String: {smiles_string}")
    print("---------------------------------------------")

if __name__ == '__main__':
    get_compound_4_details()