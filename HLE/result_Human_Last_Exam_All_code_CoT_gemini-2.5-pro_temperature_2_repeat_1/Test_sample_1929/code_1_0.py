import sys

# The user is asked to run the following python code to solve the problem.
# First, you need to install rdkit library if you have not.
# You can install it by running 'pip install rdkit' in your shell.

def get_diels_alder_product():
    """
    Determines and displays the product of the Diels-Alder reaction between
    butadiene and 1,1-dichloro-2,2-difluoroethene.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    except ImportError:
        print("Error: The 'rdkit' library is required to run this script.")
        print("Please install it by running 'pip install rdkit' in your terminal.")
        sys.exit(1)

    # Step 1: Define reactants and the major product by their SMILES strings.
    # SMILES (Simplified Molecular Input Line Entry System) is a textual
    # representation of chemical structures.

    # Reactant 1: Butadiene
    butadiene_smiles = 'C=CC=C'
    
    # Reactant 2: 1,1-dichloro-2,2-difluoroethene
    dienophile_smiles = 'C(F)(F)=C(Cl)Cl'
    
    # The major product is 4,4-dichloro-5,5-difluorocyclohexene.
    # This regioselectivity is predicted by Frontier Molecular Orbital (FMO) theory.
    product_smiles = 'C1=CCC(Cl)(Cl)C1(F)F'

    # Step 2: Create RDKit molecule objects from the SMILES strings.
    mol_butadiene = Chem.MolFromSmiles(butadiene_smiles)
    mol_dienophile = Chem.MolFromSmiles(dienophile_smiles)
    mol_product = Chem.MolFromSmiles(product_smiles)

    # A check to ensure molecules were created correctly.
    if not all([mol_butadiene, mol_dienophile, mol_product]):
        print("Error: Could not parse one or more of the molecular SMILES strings.")
        sys.exit(1)

    # Step 3: Calculate the chemical formula for each molecule.
    formula_butadiene = CalcMolFormula(mol_butadiene)
    formula_dienophile = CalcMolFormula(mol_dienophile)
    formula_product = CalcMolFormula(mol_product)

    # Step 4: Display the reaction details and product information.
    print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.")
    print("\nThe balanced chemical equation is:")
    
    # Printing with stoichiometric coefficients to show "each number in the final equation"
    print(f"1 {formula_butadiene} + 1 {formula_dienophile} -> 1 {formula_product}")
    
    print("\n--- Product Information ---")
    print(f"Name: 4,4-dichloro-5,5-difluorocyclohexene")
    print(f"Formula: {formula_product}")
    print(f"SMILES: {product_smiles}")


if __name__ == '__main__':
    get_diels_alder_product()
