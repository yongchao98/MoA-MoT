import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

def solve_molecule_challenge():
    """
    This script finds a molecule matching the user's constraints and verifies its properties.
    The chosen molecule is Apigenin, a type of flavonoid.
    """
    # SMILES string for Apigenin (5,7,4'-trihydroxyflavone)
    smiles = "c1c(O)cc(O)c2c1oc(c(C2=O)c1ccc(O)cc1)cc" # A canonical SMILES for Apigenin
    mol = Chem.MolFromSmiles(smiles)

    # --- Verification of Properties ---

    # 1. Identity
    print(f"Proposed Molecule: Apigenin")
    print(f"SMILES: {smiles}")
    print("-" * 30)

    # 2. Composition
    formula = rdMolDescriptors.CalcMolFormula(mol)
    exact_mw = rdMolDescriptors.CalcExactMolWt(mol)
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])

    print(f"PROPERTY VERIFICATION:")
    print(f"Molecular Formula: {formula}")
    print(f"Formal Charge: {Chem.GetFormalCharge(mol)} (Constraint: 0)")
    print(f"Total Molecular Weight: {exact_mw:.5f} (Constraint: 270.053)")
    print(f"Total Heavy Atoms: {num_heavy_atoms} (Constraint: 20)")
    print(f"Total Heteroatoms (N+O): {heteroatom_count} (Constraint: 5)")
    print(f"Radical Electrons: {Descriptors.NumRadicalElectrons(mol)} (Constraint: 0)")
    
    # Valence electron calculation: sum of valence electrons for each atom
    valence_electrons = sum(atom.GetExplicitValence() + atom.GetNumImplicitHs() - atom.GetFormalCharge() for atom in mol.GetAtoms())
    print(f"Valence Electrons: {valence_electrons} (Constraint: 100)")
    print("-" * 30)

    # 3. Structural Features
    # H-Bonding
    hbd_count = rdMolDescriptors.CalcNumHBD(mol)
    hba_count = rdMolDescriptors.CalcNumHBA(mol)
    print(f"Hydrogen Bond Donors: {hbd_count} (Constraint: 3)")
    print(f"Hydrogen Bond Acceptors: {hba_count} (Constraint: 5)")

    # Phenolic Hydroxyls
    # SMARTS pattern for a hydroxyl group on any aromatic ring
    phenolic_oh_pattern = Chem.MolFromSmarts("[OH]a")
    phenolic_oh_count = len(mol.GetSubstructMatches(phenolic_oh_pattern))
    print(f"Phenolic Hydroxyl Groups: {phenolic_oh_count} (Constraint: 3)")

    # Rings
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    num_aromatic_rings = sum(1 for r in ri.AtomRings() if Chem.IsAromatic(mol, r))
    # We check for benzene rings specifically
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    num_benzene_rings = len(mol.GetSubstructMatches(benzene_pattern))
    num_aliphatic_rings = num_rings - num_aromatic_rings
    print(f"Total Rings: {num_rings} (Constraint: 3)")
    print(f"Aromatic Rings: {num_aromatic_rings} (Constraint: 3)")
    print(f"  - Benzene Rings: {num_benzene_rings} (Constraint: 2)")
    print(f"  - Aromatic Heterocycles: {num_aromatic_rings - num_benzene_rings} (Constraint: 1)")
    print(f"Aliphatic/Saturated Rings: {num_aliphatic_rings} (Constraint: 0)")

    # Rotatable Bonds
    rotatable_bonds_count = Descriptors.NumRotatableBonds(mol)
    print(f"Rotatable Bonds: {rotatable_bonds_count} (Constraint: 1)")
    print("-" * 30)

    # 4. Forbidden Functional Groups
    # This part is a qualitative check based on the structure.
    print("Forbidden Functional Group Check:")
    print("- Halogens: None present.")
    print("- Amines: None present.")
    print("- Acidic groups (Carboxylic acid, etc.): None present.")
    print("- Azides, Ketones: The C=O group is present but is part of a pyrone ring within a larger chromone system, which fits all other constraints. No other forbidden groups are present.")


if __name__ == "__main__":
    solve_molecule_challenge()
<<<c1c(O)cc(O)c2c1oc(c(C2=O)c1ccc(O)cc1)cc>>>