# To run this code, you first need to install the RDKit library.
# You can do this by running: pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem

def predict_ortho_methylation_product():
    """
    This function predicts the product of a directed ortho-metalation
    and subsequent methylation reaction on N,N-diethyl-3-dimethylaminobenzamide.
    """
    # Step 1: Define the reactant molecule using its SMILES string.
    # The reactant is N,N-diethyl-3-dimethylaminobenzamide.
    reactant_smiles = 'CN(C)c1cccc(C(=O)N(CC)CC)c1'
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)

    print("--- Chemical Reaction ---")
    print("Reactant: N,N-diethyl-3-dimethylaminobenzamide")
    print("Reagents: 1) sec-BuLi, TMEDA, THF  2) CH3I")
    print("-------------------------")

    # Step 2: Define the chemical transformation using a reaction SMARTS string.
    # The reaction targets the most acidic proton, which is on the carbon (C2)
    # situated between the strong amide directing group (at C1) and the weaker
    # amine directing group (at C3).
    # The SMARTS string finds this specific C-H bond and replaces the H with a methyl group (C).
    # [c:1](C(=O)N)-[cH1:2]-[c:3](N(C)C) >> [c:1](C(=O)N)-[c:2](C)-[c:3](N(C)C)
    reaction_smarts = '[c:1](C(=O)N)-[cH1:2]-[c:3](N(C)C)>>[c:1](C(=O)N)-[c:2](C)-[c:3](N(C)C)'
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # Step 3: Run the reaction on the reactant molecule.
    # The result is a tuple containing tuples of product molecules.
    products = rxn.RunReactants((reactant_mol,))

    # Step 4: Extract and display the primary product.
    if products and products[0]:
        product_mol = products[0][0]
        # Clean up the molecule's representation
        Chem.SanitizeMol(product_mol)
        
        # Convert the product molecule back to a SMILES string
        product_smiles = Chem.MolToSmiles(product_mol)
        
        # The known IUPAC name of the product
        product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"

        print("\n--- Predicted Product ---")
        print("Product Name:   ", product_name)
        print("Product SMILES: ", product_smiles)

    else:
        print("Reaction model failed to produce a product.")

if __name__ == '__main__':
    predict_ortho_methylation_product()