# The user needs to install the rdkit and pubchempy libraries first.
# You can install them by running this command in your terminal or command prompt:
# pip install rdkit pubchempy

from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
import sys

def find_product_iupac_name():
    """
    This function determines the product of the specified chemical reaction and prints its IUPAC name.
    
    The reaction is between 1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.
    
    Chemical Rationale:
    1.  The Grignard reagent (phenyl magnesium bromide) is a potent nucleophile.
    2.  The substrate has three halogen leaving groups (I, Br, Br).
    3.  Reactivity of aryl-halogen bonds towards Grignard reagents is I >> Br.
    4.  Therefore, a selective substitution of the iodine atom by a phenyl group will occur.
    5.  The less reactive C-Br bonds will not react under these (uncatalyzed) conditions.
    6.  The final product is 1,3-dibromo-2-phenylbenzene.
    
    This script models this chemical transformation and uses an external database (PubChem) to find the correct IUPAC name.
    """
    # SMILES representation of the reactant 1,3-dibromo-2-iodobenzene
    reactant_smiles = "c1(Br)c(I)c(Br)ccc1"
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)

    if not reactant_mol:
        print("Error: Invalid reactant SMILES string.")
        return

    # The reaction is defined using SMARTS: an aromatic carbon bonded to iodine ([c:1]-I)
    # is transformed into the same carbon bonded to a phenyl group ([c:1]-c1ccccc1).
    reaction_smarts = "[c:1]-I>>[c:1]-c1ccccc1"
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # Execute the reaction
    products = rxn.RunReactants((reactant_mol,))

    if not products or not products[0]:
        print("Error: The reaction did not yield a product.")
        return
        
    # Get the first product molecule from the result
    product_mol = products[0][0]
    Chem.SanitizeMol(product_mol) # Clean up the molecule's structure
    
    # Generate a canonical SMILES string for the product to ensure a standard representation
    product_smiles = Chem.MolToSmiles(product_mol, canonical=True)
    
    print(f"Reactant: 1,3-dibromo-2-iodobenzene")
    print(f"Reaction: Substitution of iodine with a phenyl group.")
    print(f"Product Structure (SMILES): {product_smiles}")
    
    try:
        # Use the SMILES string to search PubChem for the compound
        print("\nSearching for IUPAC name in PubChem database...")
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        if compounds:
            iupac_name = compounds[0].iupac_name
            print("\n--- Final Answer ---")
            print(f"The IUPAC name of the product is: {iupac_name}")
            
            # As requested, outputting the numbers from the IUPAC name
            numbers = [char for char in iupac_name if char.isdigit()]
            print(f"The numbers in the name are: {', '.join(numbers)}")

        else:
            print("Could not find the compound in PubChem. The predicted name is 1,3-dibromo-2-phenylbenzene.")

    except Exception as e:
        print(f"\nAn error occurred while communicating with PubChem: {e}", file=sys.stderr)
        print("Based on chemical principles, the predicted IUPAC name is: 1,3-dibromo-2-phenylbenzene")

if __name__ == "__main__":
    find_product_iupac_name()