# To run this code, you need to have the RDKit library installed.
# You can install it using pip: pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem

def get_molecular_formula(mol):
    """Calculates the molecular formula from an RDKit molecule object."""
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def solve_reaction():
    """
    Solves the reaction between butadiene and 1,1-dichloro-2,2-difluoroethene.
    """
    # Step 1: Define reactants using SMILES strings
    # Butadiene: CH2=CH-CH=CH2
    diene_smiles = 'C=CC=C'
    # 1,1-dichloro-2,2-difluoroethene: CCl2=CF2
    dienophile_smiles = 'C(Cl)(Cl)=C(F)F'

    diene_mol = Chem.MolFromSmiles(diene_smiles)
    dienophile_mol = Chem.MolFromSmiles(dienophile_smiles)

    if not diene_mol or not dienophile_mol:
        print("Error: Invalid SMILES string for reactants.")
        return

    # Step 2: Define the Diels-Alder reaction using reaction SMARTS
    # The pattern is [diene].[dienophile]>>[product]
    # The numbers (:1, :2, etc.) map atoms from reactants to products.
    reaction_smarts = '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1'
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # Step 3: Run the reaction
    reactants = (diene_mol, dienophile_mol)
    products_sets = rxn.RunReactants(reactants)

    if not products_sets:
        print("The reaction did not yield any products with the given SMARTS.")
        return

    # Step 4: Get the product and its information
    product_mol = products_sets[0][0]  # The main product is in the first set
    Chem.SanitizeMol(product_mol)

    # Manually determined IUPAC name based on chemical principles
    product_iupac_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"
    
    # Get molecular formulas for the equation
    diene_formula = get_molecular_formula(diene_mol)
    dienophile_formula = get_molecular_formula(dienophile_mol)
    product_formula = get_molecular_formula(product_mol)

    # Step 5: Print the results as a final equation
    print("The reaction is a Diels-Alder cycloaddition.")
    print("The final balanced chemical equation is:")
    # The coefficient for each molecule is 1
    print(f"1 {diene_formula} + 1 {dienophile_formula} -> 1 {product_formula}")
    print("\nThe product is:")
    print(product_iupac_name)

    # Step 6: Output each number from the product name, as requested
    print("\nThe numbers in the product's IUPAC name are:")
    print("4")
    print("4")
    print("5")
    print("5")
    print("1")

solve_reaction()