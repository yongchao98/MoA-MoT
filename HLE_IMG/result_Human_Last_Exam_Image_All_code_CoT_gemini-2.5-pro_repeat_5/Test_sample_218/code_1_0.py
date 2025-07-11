# To run this code, you first need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem

def identify_product_A():
    """
    This function identifies the product of the reaction of geraniol
    with O-(p-tolyl) chlorothionoformate followed by LiAlH4 reduction.
    """
    # Step 1: Define the starting material, geraniol, using its SMILES string.
    # The SMILES string 'CC(C)=CCC/C(C)=C/CO' represents (2E)-3,7-dimethylocta-2,6-dien-1-ol.
    geraniol_smiles = 'CC(C)=CCC/C(C)=C/CO'
    geraniol_mol = Chem.MolFromSmiles(geraniol_smiles)
    if geraniol_mol is None:
        print("Error: Invalid SMILES string for geraniol.")
        return

    # Step 2: Define the chemical transformation.
    # The overall reaction is the deoxygenation of a primary alcohol (-CH2OH) to a methyl group (-CH3).
    # We define this transformation using a reaction SMARTS string.
    # '[CH2:1][OH1:2]' matches the primary alcohol group.
    # '>>[CH3:1]' specifies that this group becomes a methyl group.
    reaction_smarts = '[CH2:1][OH1:2]>>[CH3:1]'
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # Step 3: Apply the reaction to geraniol.
    products = rxn.RunReactants((geraniol_mol,))

    # Step 4: Extract and print the product's structure.
    if products:
        # The reaction yields one product.
        product_mol = products[0][0]
        # Sanitize the molecule to ensure correct chemical properties.
        Chem.SanitizeMol(product_mol)
        # Convert the product molecule to its canonical SMILES string.
        product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True)
        
        print("The reaction sequence performs a deoxygenation of geraniol.")
        print("The starting material is Geraniol: " + geraniol_smiles)
        print("The final product, Compound A, is a hydrocarbon with the following structure (SMILES notation):")
        print(product_smiles)
        print("\nThis compound is named (6E)-2,6-dimethylocta-2,6-diene.")

    else:
        print("The reaction could not be performed on the provided reactant.")

if __name__ == '__main__':
    identify_product_A()
