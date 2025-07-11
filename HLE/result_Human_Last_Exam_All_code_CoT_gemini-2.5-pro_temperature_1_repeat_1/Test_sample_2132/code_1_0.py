# First, ensure you have RDKit installed:
# pip install rdkit-pypi

from rdkit import Chem

def calculate_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: The product of the Favorskii rearrangement of 2-chlorocyclohexanone
    # is cyclopentanecarboxylic acid. We represent it with a SMILES string.
    smiles_product = "O=C(O)C1CCCC1"
    molecule = Chem.MolFromSmiles(smiles_product)

    if molecule is None:
        print("Error: Could not create molecule from SMILES string.")
        return

    # Step 2: Get the number of non-hydrogen atoms and bonds for the complexity formula.
    # RDKit's GetNumAtoms() and GetNumBonds() on a molecule from SMILES
    # give the counts for the hydrogen-suppressed graph.
    num_atoms = molecule.GetNumAtoms()
    num_bonds = molecule.GetNumBonds()

    # Step 3: Calculate the complexity using the formula: (N_atoms * N_bonds) / 2
    complexity = (num_atoms * num_bonds) / 2

    # Step 4: Print the details of the calculation as requested.
    print("Molecule: Cyclopentanecarboxylic Acid")
    print("Complexity Formula: (Number of non-hydrogen atoms * Number of non-hydrogen bonds) / 2")
    print("\n--- Calculation ---")
    print(f"Number of non-hydrogen atoms (N_atoms) = {num_atoms}")
    print(f"Number of non-hydrogen bonds (N_bonds) = {num_bonds}")
    print(f"Complexity = ({num_atoms} * {num_bonds}) / 2")
    print(f"Final Böttcher Molecular Complexity = {complexity}")

if __name__ == "__main__":
    calculate_complexity()