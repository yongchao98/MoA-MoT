# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

def solve_reaction():
    """
    This function identifies Compound A from the described reaction,
    calculates its properties, and prints the information.
    """
    # SMILES string for the final product, Compound A: 2-(pyridin-2-yl)-1H-isoindole-3-carbonitrile
    # This structure is derived from the Joullié-Ugi three-component reaction mechanism.
    product_smiles = "n1ccccc1N1C=C(C#N)c2c1cccc2"
    
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(product_smiles)

    if mol:
        # Calculate molecular formula and molecular weight
        mol_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        mol_weight = Descriptors.MolWt(mol)

        print("--- Analysis of Compound A ---")
        print(f"Product Name: 2-(pyridin-2-yl)-1H-isoindole-3-carbonitrile")
        print(f"Molecular Formula: {mol_formula}")
        print(f"Molecular Weight: {mol_weight:.4f} g/mol")
        print("\n--- Net Reaction Equation ---")
        
        # The equation shows the main organic components and stoichiometry
        # Reactants: 1 (2-aminopyridine) + 1 (o-phthalaldehyde) + 1 (HCN source)
        # Products: 1 (Compound A) + 2 (Water)
        equation = "1 C5H6N2 + 1 C8H6O2 + 1 HCN -> 1 C14H9N3 + 2 H2O"
        print(equation)
        print("\nNote: The equation is balanced for atoms contributing to the final structure, with HCN representing the cyanide source and water as the byproduct.")

    else:
        print("Error: Could not generate molecule from SMILES string.")

# Execute the function
solve_reaction()
