# Import necessary libraries from the RDKit package
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def design_and_verify_molecule():
    """
    Designs a molecule based on a specific set of constraints and verifies its properties.
    The strategy is to build a complex cage-like molecule by programmatically adding atoms and bonds.
    """

    # 1. Initialize an editable molecule object
    mol = Chem.RWMol()

    # 2. Add the atoms: 6 Oxygens and 12 Carbons, to match the C12O6 heavy-atom skeleton.
    for _ in range(6): mol.AddAtom(Chem.Atom(8))  # Add 6 Oxygen atoms
    for _ in range(12): mol.AddAtom(Chem.Atom(6)) # Add 12 Carbon atoms

    # Atom indices for clarity in bonding:
    # Oxygens: 0, 1, 2, 3, 4, 5
    # Carbons: 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
    
    # 3. Add bonds to create a specific tricyclic ketone structure
    # This structure is a modification of an 18-crown-6 ether, made rigid by cross-links.
    
    # Main ring bonds (ether and C-C linkages)
    mol.AddBond(6, 1, Chem.BondType.SINGLE)       # C-O
    mol.AddBond(1, 8, Chem.BondType.SINGLE)       # O-C
    mol.AddBond(8, 9, Chem.BondType.SINGLE)       # C-C
    mol.AddBond(9, 2, Chem.BondType.SINGLE)       # C-O
    mol.AddBond(2, 10, Chem.BondType.SINGLE)      # O-C
    mol.AddBond(10, 11, Chem.BondType.SINGLE)     # C-C
    mol.AddBond(11, 3, Chem.BondType.SINGLE)      # C-O
    mol.AddBond(3, 13, Chem.BondType.SINGLE)      # O-C
    mol.AddBond(13, 12, Chem.BondType.SINGLE)     # C-C
    mol.AddBond(12, 4, Chem.BondType.SINGLE)      # C-O
    mol.AddBond(4, 15, Chem.BondType.SINGLE)      # O-C
    mol.AddBond(15, 14, Chem.BondType.SINGLE)     # C-C
    mol.AddBond(14, 5, Chem.BondType.SINGLE)      # C-O
    mol.AddBond(5, 17, Chem.BondType.SINGLE)      # O-C
    mol.AddBond(17, 16, Chem.BondType.SINGLE)     # C-C
    mol.AddBond(16, 0, Chem.BondType.SINGLE)      # C-O

    # Bond for the carbonyl group (C=O)
    mol.AddBond(6, 0, Chem.BondType.DOUBLE)      # C=O
    
    # Bonds to create the rigid, tricyclic cage structure (cross-links)
    mol.AddBond(7, 16, Chem.BondType.SINGLE)      # C-C cross-link 1
    mol.AddBond(7, 11, Chem.BondType.SINGLE)      # C-C cross-link 2
    
    # 4. Finalize the molecule from its editable form
    final_mol = mol.GetMol()
    
    # Clean up the molecule's representation and add hydrogens
    Chem.SanitizeMol(final_mol)
    final_mol = Chem.AddHs(final_mol)
    
    # 5. Generate and print the canonical SMILES string
    smiles = Chem.MolToSmiles(final_mol, isomericSmiles=False)
    
    # --- Verification Step ---
    # The following code verifies the molecule against the provided constraints.
    # The constraints are listed with the calculated value for the generated molecule.
    
    print("--- Verification of Molecular Properties ---")
    
    # Molecular Formula and Weight
    formula = rdMolDescriptors.CalcMolFormula(final_mol)
    mw = Descriptors.ExactMolWt(final_mol)
    print(f"Molecular Formula: {formula} (Expected: C12H18O6)")
    print(f"Molecular Weight: {mw:.2f} g/mol (Target: 258.11)")
    
    # Atom Counts
    heavy_atoms = Descriptors.HeavyAtomCount(final_mol)
    heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(final_mol)
    print(f"Total Heavy Atoms: {heavy_atoms} (Target: 18)")
    print(f"Total Heteroatoms: {heteroatoms} (Target: 6)")
    
    # Electron Count
    valence_electrons = Descriptors.NumValenceElectrons(final_mol)
    print(f"Valence Electrons: {valence_electrons} (Target: 102)")

    # Ring Information
    ring_info = final_mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    print(f"Total Rings: {num_rings} (Target: 3)")

    # H-Bond Donors/Acceptors
    hbd = Descriptors.NumHDonors(final_mol)
    hba = Descriptors.NumHAcceptors(final_mol)
    print(f"Hydrogen Bond Donors: {hbd} (Target: 0)")
    print(f"Hydrogen Bond Acceptors: {hba} (Target: 6)")
    
    # Structural Features
    rot_bonds = Descriptors.NumRotatableBonds(final_mol)
    is_aromatic = any(b.GetIsAromatic() for b in final_mol.GetBonds())
    print(f"Rotatable Bonds: {rot_bonds} (Target: 0)")
    print(f"Contains Aromatic Rings: {is_aromatic} (Target: False)")
    
    # Functional Group Checks (SMARTS patterns for key groups)
    num_carbonyls = len(final_mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[O]')))
    num_esters = len(final_mol.GetSubstructMatches(Chem.MolFromSmarts('[#6][#6](=O)O[#6]')))
    print(f"Number of Carbonyls: {num_carbonyls} (Target: >=1)")
    print(f"Number of Esters: {num_esters} (Target: 0)")

    # The problem demands the final answer in a specific format.
    # The SMILES string is the solution to the molecular design task.
    print("\nFinal molecule represented in SMILES format:")
    print(f"<<<{smiles}>>>")

if __name__ == "__main__":
    design_and_verify_molecule()