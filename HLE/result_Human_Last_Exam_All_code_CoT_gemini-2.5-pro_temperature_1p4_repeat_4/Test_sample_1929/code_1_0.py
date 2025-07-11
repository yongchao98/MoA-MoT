# The user needs to have rdkit installed. This can be done by running:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def solve_reaction():
    """
    Solves the Diels-Alder reaction between butadiene and 
    1,1-dichloro-2,2-difluoroethene.
    """
    # 1. Define reactants using SMILES strings
    # Butadiene: CH2=CH-CH=CH2
    diene_smiles = "C=CC=C"
    diene_name = "1,3-Butadiene"

    # 1,1-dichloro-2,2-difluoroethene: CCl2=CF2
    dienophile_smiles = "C(Cl)(Cl)=C(F)F"
    dienophile_name = "1,1-dichloro-2,2-difluoroethene"

    # Create RDKit molecule objects
    diene_mol = Chem.MolFromSmiles(diene_smiles)
    dienophile_mol = Chem.MolFromSmiles(dienophile_smiles)

    if not diene_mol or not dienophile_mol:
        print("Error: Could not parse reactant SMILES.")
        return

    # 2. Define the Diels-Alder reaction using a reaction SMARTS string
    # [diene:1].[dienophile:2]>>[product:1]
    # Maps atoms from reactants to products to form the cyclohexene ring
    rxn_smarts = "[c:1]=[c:2][c:3]=[c:4].[c:5]=[c:6]>>[c:1]1[c:2]=[c:3][c:4][c:6][c:5]1"
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    # 3. Run the reaction
    reactants = (diene_mol, dienophile_mol)
    products = rxn.RunReactants(reactants)

    # 4. Analyze and print the result
    print("The reaction is a Diels-Alder [4+2] cycloaddition.")
    print("-" * 50)
    
    if products:
        # The reaction yields one major product
        product_mol = products[0][0]
        Chem.SanitizeMol(product_mol)

        # Determine product properties
        product_name = "3,3-dichloro-4,4-difluorocyclohexene"
        product_formula = CalcMolFormula(product_mol)
        
        # Print the final equation with coefficients and names
        print("Final Equation:")
        print(f"1 * ({diene_name}) + 1 * ({dienophile_name}) -> 1 * ({product_name})")
        print("\nProduct Details:")
        print(f"  Name: {product_name}")
        print(f"  Molecular Formula: {product_formula}")

    else:
        print("Reaction simulation failed to produce a product.")

# Execute the function
solve_reaction()
<<<3,3-dichloro-4,4-difluorocyclohexene>>>