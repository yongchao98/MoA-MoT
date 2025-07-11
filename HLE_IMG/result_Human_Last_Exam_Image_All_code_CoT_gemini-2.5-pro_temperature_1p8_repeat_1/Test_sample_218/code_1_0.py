# The script requires the rdkit-pypi package.
# You can install it using pip: pip install rdkit-pypi

def solve_chemistry_problem():
    """
    This script identifies Compound A and calculates its chemical properties.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    except ImportError:
        print("RDKit library not found. Please install it using 'pip install rdkit-pypi' to run this script.")
        print("\nBased on chemical principles, Compound A is (2E)-3,7-dimethylocta-2,6-diene.")
        return

    # --- Step-by-step Chemical Analysis ---
    print("### Analysis of the Reaction ###")
    print("1. The starting material is geraniol, a primary allylic alcohol.")
    print("2. In the first step, the alcohol (-OH) is converted to an O-aryl thionocarbonate.")
    print("3. In the second step, reduction with LiAlH4 removes the thionocarbonate group and replaces it with a hydrogen atom.")
    print("4. This is a deoxygenation reaction, transforming the -CH2OH group into a -CH3 group.")
    
    # --- Identification of Compound A ---
    # The SMILES string for (2E)-3,7-dimethylocta-2,6-diene is used.
    # Structure: (CH3)2C=CH-CH2-CH2-C(CH3)=CH-CH3
    product_smiles = "CC(C)=CCC/C(C)=C/C"
    
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(product_smiles)

    if mol is None:
        print("Error: Could not generate molecule from SMILES string.")
        return

    # --- Properties of Compound A ---
    print("\n### Properties of Compound A ###")

    # Name, Formula, and Mass
    iupac_name = "(2E)-3,7-dimethylocta-2,6-diene"
    molecular_formula = CalcMolFormula(mol)
    average_mass = Descriptors.MolWt(mol)

    print(f"IUPAC Name: {iupac_name}")
    print(f"SMILES: {product_smiles}")
    print(f"Molecular Formula: {molecular_formula}")
    
    # --- Molar Mass Calculation ---
    # Atomic weights (average)
    atomic_mass_C = 12.011
    atomic_mass_H = 1.008
    
    # Number of atoms from formula C10H18
    num_C = 10
    num_H = 18

    print("\nMolar Mass Calculation:")
    print(f"The calculation is based on the number of Carbon and Hydrogen atoms in the formula {molecular_formula}.")
    print(f"Equation: ({num_C} * mass_of_C) + ({num_H} * mass_of_H)")
    print(f"Result: ({num_C} * {atomic_mass_C}) + ({num_H} * {atomic_mass_H}) = {average_mass:.3f} g/mol")

if __name__ == '__main__':
    solve_chemistry_problem()
