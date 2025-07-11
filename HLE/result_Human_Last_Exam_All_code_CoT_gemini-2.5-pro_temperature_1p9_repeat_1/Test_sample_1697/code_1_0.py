# RDKit is a cheminformatics toolkit. If you don't have it, you can install it:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

def get_product_info():
    """
    This function determines the product of the described reaction
    and prints its information.
    """
    # Step 1: Reasoning (as explained above)
    # The reaction is a directed ortho-metalation at C2, followed by methylation.
    # Starting Material: N,N-diethyl-3-dimethylaminobenzamide
    # SMILES: CCN(CC)C(=O)c1cccc(c1)N(C)C
    # Product: N,N-diethyl-2-methyl-3-dimethylaminobenzamide

    # Step 2: Define the product using its SMILES string
    product_smiles = "CCN(CC)C(=O)c1c(C)c(N(C)C)ccc1"
    
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(product_smiles)

    if mol is None:
        print("Error: Could not generate molecule from SMILES string.")
        return

    # To get a reasonable IUPAC-like name, we can try to generate it,
    # but for complex molecules, manual naming is more reliable.
    # The systematic name is derived from the chemical reasoning above.
    product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"

    # Calculate molecular formula, including all numbers for atoms
    molecular_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    
    print("Final Product Information:")
    print(f"IUPAC Name: {product_name}")
    print(f"SMILES String: {product_smiles}")
    print(f"Molecular Formula: {molecular_formula}")

# Execute the function to get the answer
get_product_info()
