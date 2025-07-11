# The 'pubchempy' library is required. If not installed, run: pip install pubchempy
import pubchempy as pcp
import sys

def find_reaction_product():
    """
    This script determines the product of a chemical reaction and fetches its data.
    Reaction: N,N-diethyl-3-dimethylaminobenzamide reacts first with sec-BuLi/TMEDA
              and then with methyl iodide.
    The predicted product is N,N-diethyl-2-methyl-3-dimethylaminobenzamide.
    """

    # The SMILES (Simplified Molecular Input Line Entry System) string for the final product.
    # Structure: A benzene ring with -CON(CC)2 at C1, -CH3 at C2, and -N(C)2 at C3.
    product_smiles = "CN(C)c1c(C)cccc1C(=O)N(CC)CC"
    
    print("Step 1: Directed Ortho-Metalation")
    print("The strong base (sec-BuLi with TMEDA) removes a proton from the position ortho (adjacent) to the directing groups.")
    print("Both the amide (-CONEt2) and amine (-NMe2) groups direct the lithiation to position 2.\n")
    
    print("Step 2: Electrophilic Quench")
    print("The resulting carbanion at position 2 attacks the methyl group from methyl iodide (CH3I).\n")
    
    print("--- Final Product ---")
    
    try:
        # Search PubChem for the compound using its SMILES string.
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        
        if not compounds:
            print("Compound not found on PubChem. Manually providing predicted information.")
            print("Predicted IUPAC Name: N,N-diethyl-2-methyl-3-(dimethylamino)benzamide")
            return

        product = compounds[0]
        
        # Output the compound details
        print(f"IUPAC Name: {product.iupac_name}")
        print(f"Molecular Formula: {product.molecular_formula}")
        print(f"Molecular Weight: {product.molecular_weight} g/mol")
        print(f"PubChem CID: {product.cid}")

    except Exception as e:
        # Handle potential connection errors to the PubChem API
        print(f"An error occurred while fetching data from PubChem: {e}", file=sys.stderr)
        print("Falling back to predicted information:")
        print("Predicted IUPAC Name: N,N-diethyl-2-methyl-3-(dimethylamino)benzamide")


# Run the function
find_reaction_product()