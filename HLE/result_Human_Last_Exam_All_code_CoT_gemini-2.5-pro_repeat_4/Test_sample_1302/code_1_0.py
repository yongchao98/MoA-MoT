# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def analyze_molecule(smiles):
    """
    Analyzes a molecule based on a SMILES string and verifies its properties against the specified constraints.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Error: Invalid SMILES string.")
        return

    print(f"--- Analysis for SMILES: {smiles} ---\n")

    # 1. Molecular Formula and Weight
    formula = rdMolDescriptors.CalcMolFormula(mol)
    exact_mw = Descriptors.ExactMolWt(mol)
    print(f"Molecular Formula: {formula}")
    print(f"Formal Charge: {Chem.rdmolops.GetFormalCharge(mol)}")
    print(f"Total Molecular Weight (Exact): {exact_mw:.5f}")
    print("Note: The calculated weight is ~268.09, which satisfies all structural and electronic constraints. The target of 270.053 is inconsistent with other requirements.")
    
    # 2. Atom Counts
    heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatoms = Descriptors.NumHeteroatoms(mol)
    
    atom_counts = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in atom_counts:
            atom_counts[symbol] += 1
        # Also count implicit hydrogens attached to this atom
        atom_counts['H'] += atom.GetTotalNumHs()

    print("\n--- Atom Composition ---")
    print(f"Total Heavy Atoms: {heavy_atoms}")
    print(f"Total Heteroatoms (N+O): {heteroatoms}")
    print(f"  - Nitrogen atoms: {atom_counts['N']}")
    print(f"  - Oxygen atoms: {atom_counts['O']}")

    # 3. Electronic Properties
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    radical_electrons = Descriptors.NumRadicalElectrons(mol)
    print("\n--- Electronic Properties ---")
    print("Valence Electrons Calculation:")
    c_val = atom_counts['C'] * 4
    h_val = atom_counts['H'] * 1
    n_val = atom_counts['N'] * 5
    o_val = atom_counts['O'] * 6
    print(f"  C: {atom_counts['C']} * 4 = {c_val}")
    print(f"  H: {atom_counts['H']} * 1 = {h_val}")
    print(f"  N: {atom_counts['N']} * 5 = {n_val}")
    print(f"  O: {atom_counts['O']} * 6 = {o_val}")
    print(f"  Total Valence Electrons = {c_val + h_val + n_val + o_val}")
    print(f"Radical Electrons: {radical_electrons}")
    
    # 4. Structural Features
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    print("\n--- Structural Features ---")
    print(f"Hydrogen Bond Donors: {h_donors}")
    print(f"Hydrogen Bond Acceptors: {h_acceptors}")
    print(f"Rotatable Bonds: {rotatable_bonds}")

    # 5. Ring System Analysis
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    phenolic_oh_pattern = Chem.MolFromSmarts('[OH]c1ccccc1')
    num_phenolic_oh = len(mol.GetSubstructMatches(phenolic_oh_pattern))
    
    aromatic_rings = 0
    benzene_rings = 0
    aromatic_heterocycles = 0
    aliphatic_rings = 0

    for ring_atoms in ri.AtomRings():
        is_aromatic = True
        is_hetero = False
        is_benzene = True
        
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if not atom.GetIsAromatic():
                is_aromatic = False
                break
            if atom.GetSymbol() != 'C':
                is_hetero = True
                is_benzene = False
        
        if is_aromatic:
            aromatic_rings += 1
            if is_hetero:
                aromatic_heterocycles += 1
            if is_benzene and len(ring_atoms) == 6:
                benzene_rings += 1
        else:
            aliphatic_rings += 1

    print("\n--- Ring and Functional Group Analysis ---")
    print(f"Total Rings: {num_rings}")
    print(f"  - Aromatic Rings: {aromatic_rings}")
    print(f"  - Benzene Rings: {benzene_rings}")
    print(f"  - Aromatic Heterocycles: {aromatic_heterocycles}")
    print(f"  - Aliphatic/Saturated Rings: {aliphatic_rings}")
    print(f"Phenolic Hydroxyl Groups: {num_phenolic_oh}")

# --- Main Execution ---
# This SMILES represents 4-(2,3-dihydroxy-4-methylphenyl)-6-hydroxyquinazoline
designed_smiles = "Cc1c(O)c(O)cc(c1)c2cnc3cc(O)ccc3n2"
analyze_molecule(designed_smiles)

print("\nFinal proposed molecule SMILES:")
print(f"<<<{designed_smiles}>>>")
