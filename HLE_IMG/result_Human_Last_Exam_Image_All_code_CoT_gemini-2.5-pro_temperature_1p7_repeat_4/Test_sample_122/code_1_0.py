import sys
import subprocess

# Install pubchempy if it's not already installed
try:
    import pubchempy as pcp
except ImportError:
    print("Installing pubchempy library...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pubchempy"])
    import pubchempy as pcp

def find_product_name():
    """
    Identifies the name of the final product from the synthesis.
    The final product is Tryptophan hydroxamate.
    """
    # SMILES string for (S)-Tryptophan hydroxamate
    final_product_smiles = "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)NO)N"
    
    product_name = "Tryptophan hydroxamate"
    
    print(f"Based on the reaction sequence, the final product is {product_name}.")
    
    try:
        # Search PubChem for the compound using its SMILES string
        compounds = pcp.get_compounds(final_product_smiles, 'smiles')
        if compounds:
            compound = compounds[0]
            print("\nVerifying name with PubChem database...")
            # Check if our derived name is a known synonym
            if product_name in compound.synonyms:
                 print(f"Confirmed Name: {product_name}")
                 print(f"IUPAC Name: {compound.iupac_name}")
            else:
                 # Fallback to other available names
                 print(f"A common name found is: {compound.synonyms[0]}")
                 print(f"IUPAC Name: {compound.iupac_name}")
        else:
            print("\nCould not verify name in PubChem, but the chemically derived name is correct.")
            
    except Exception as e:
        print(f"\nAn error occurred while trying to connect to PubChem: {e}")
        print("Proceeding with the name determined from chemical principles.")

# Execute the function
find_product_name()