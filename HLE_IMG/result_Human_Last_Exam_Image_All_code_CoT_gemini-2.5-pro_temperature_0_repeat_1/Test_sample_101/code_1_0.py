# First, ensure you have the rdkit library installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def analyze_compound_A():
    """
    This function analyzes the final product (Compound A) of the described reaction.
    The reaction is a two-step synthesis:
    1. Imine formation: 3-hydroxy-pyridine-2-carbaldehyde + aniline -> imine intermediate
    2. Cyanide addition (Strecker-type reaction): imine + NaCN -> Compound A (an alpha-aminonitrile)
    """
    # The SMILES (Simplified Molecular Input Line Entry System) string for Compound A,
    # which is 2-(cyano(phenylamino)methyl)pyridin-3-ol.
    smiles_A = "N#CC(Nc1ccccc1)c2ncccc2O"

    # Create an RDKit molecule object from the SMILES string
    mol_A = Chem.MolFromSmiles(smiles_A)

    if mol_A:
        # Add explicit hydrogen atoms to the molecule graph for accurate calculations
        mol_A_with_hs = Chem.AddHs(mol_A)

        # Calculate the molecular formula
        molecular_formula = rdMolDescriptors.CalcMolFormula(mol_A_with_hs)

        # Calculate the exact molecular weight
        exact_mw = Descriptors.ExactMolWt(mol_A_with_hs)

        print("--- Analysis of Final Product: Compound A ---")
        print("The final product, Compound A, is 2-(cyano(phenylamino)methyl)pyridin-3-ol.")
        print(f"SMILES Representation: {smiles_A}")
        print(f"Molecular Formula: {molecular_formula}")
        print(f"Exact Molecular Weight: {exact_mw:.4f} g/mol")
        print("\n--- Atom Count ('Numbers in the Final Equation') ---")
        print("The molecular formula represents the number of atoms of each element in one molecule of Compound A:")
        
        # Extract and print the count for each element from the formula
        atom_counts = rdMolDescriptors.GetMolFormula(mol_A_with_hs, separateCounts=True)
        print(f"Number of Carbon (C) atoms: {atom_counts[0]}")
        print(f"Number of Hydrogen (H) atoms: {atom_counts[1]}")
        print(f"Number of Nitrogen (N) atoms: {atom_counts[3]}")
        print(f"Number of Oxygen (O) atoms: {atom_counts[2]}")

    else:
        print("Error: Could not create molecule from the provided SMILES string.")

# Execute the analysis
analyze_compound_A()