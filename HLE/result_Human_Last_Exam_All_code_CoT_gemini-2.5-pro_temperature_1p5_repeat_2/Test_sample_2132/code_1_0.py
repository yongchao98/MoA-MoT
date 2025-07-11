import math
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def solve_molecular_complexity():
    """
    Calculates the Böttcher (Bertz) Molecular Complexity for the product of
    the Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: Identify the product molecule.
    # The rearrangement product is cyclopentanecarboxylic acid.
    # We define it using its SMILES string.
    product_smiles = "C1CCC(C1)C(=O)O"
    
    print("--- Calculating Böttcher (Bertz) Molecular Complexity ---")
    print(f"Reactant: 2-chlorocyclohexanone")
    print(f"Reaction: Favorskii Rearrangement")
    print(f"Product: Cyclopentanecarboxylic acid (SMILES: {product_smiles})\n")

    # Step 2: Prepare the molecule for calculation using RDKit.
    # The calculation requires explicit hydrogen atoms.
    mol = Chem.MolFromSmiles(product_smiles)
    mol_with_hs = Chem.AddHs(mol)

    # Step 3: Explain the Bertz Complexity formula and its components.
    # The final equation is C = C_atoms + C_bonds.
    # This calculation depends on the number of atoms and bonds and their distribution.
    
    # Get the primary numbers that go into the final equation.
    num_atoms = mol_with_hs.GetNumAtoms()
    num_bonds = mol_with_hs.GetNumBonds()

    print("The Bertz Complexity formula is C = (Atom Complexity Term) + (Bond Complexity Term).")
    print("The numbers that go into this final equation are:")
    print(f"1. Total number of atoms (Na): {num_atoms}")
    print(f"2. Total number of bonds (Nb): {num_bonds}")
    print("3. Counts of topologically unique atom and bond types (determined by the library).\n")

    # Step 4: Calculate the final complexity value.
    bertz_complexity = GraphDescriptors.BertzCT(mol_with_hs)

    # Step 5: Output the final result.
    print("The final calculated Böttcher (Bertz) Molecular Complexity is:")
    print(f"Result = {bertz_complexity}")

if __name__ == "__main__":
    solve_molecular_complexity()
